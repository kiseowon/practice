"""Four input conditions across MLP references, pure cGAN, and epsilon-DDPM."""
import argparse
import copy
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace
import uuid
import joblib
import numpy as np
import pandas as pd
import torch
from torch import nn
import torch.nn.functional as F
from threadpoolctl import threadpool_limits
import first_structure_mlp as first
import pathway_structure_mlp as pathway
from pathway_pattern_mlp import PatternMLP


MEMBERSHIP=None
class Context(nn.Module):
    def __init__(self,dims,use_pathway):
        super().__init__();gene,path,cond,structure,out=dims
        template=PatternMLP(gene,cond,structure,out,MEMBERSHIP,use_pathway)
        self.gene_dim=gene
        self.expression_encoder=template.expression_encoder
        self.structure_encoder=template.structure_encoder
        self.pathway_encoder=template.pathway_encoder
        self.fusion=template.decoder[:3]
    def forward(self,source,condition,structure):
        genes=source[:,:self.gene_dim];scores=source[:,self.gene_dim:]
        return self.fusion(torch.cat([self.expression_encoder(genes),self.structure_encoder(structure),self.pathway_encoder(genes,scores),condition],1))


class Generator(nn.Module):
    def __init__(self,dims,use_pathway):
        super().__init__();self.context=Context(dims,use_pathway)
        self.head=nn.Sequential(nn.Linear(64+16,64),nn.SiLU(),nn.Linear(64,dims[-1]))
    def decode(self,context,noise):return self.head(torch.cat([context,noise],1))
    def forward(self,*args):return self.decode(self.context(*args[:3]),args[3])


class Discriminator(nn.Module):
    def __init__(self,dims,use_pathway):
        super().__init__();self.context=Context(dims,use_pathway)
        self.head=nn.Sequential(nn.Linear(64+dims[-1],128),nn.LeakyReLU(.2),nn.Linear(128,64),nn.LeakyReLU(.2),nn.Linear(64,1))
    def forward(self,source,condition,structure,expression):return self.head(torch.cat([self.context(source,condition,structure),expression],1)).squeeze(1)


class EpsilonDenoiser(nn.Module):
    def __init__(self,dims,use_pathway,steps):
        super().__init__();self.context=Context(dims,use_pathway);self.time=nn.Embedding(steps,16);self.output_dim=dims[-1]
        self.correction=nn.Sequential(nn.Linear(64+self.output_dim+16,128),nn.SiLU(),nn.Linear(128,128),nn.SiLU(),nn.Linear(128,self.output_dim))
        nn.init.zeros_(self.correction[-1].weight);nn.init.zeros_(self.correction[-1].bias)
    def decode(self,context,noisy,t):
        # A neural residual connection, not a separate prediction/post-correction model.
        return noisy+self.correction(torch.cat([context,noisy,self.time(t)],1))
    def forward(self,source,condition,structure,noisy,t):return self.decode(self.context(source,condition,structure),noisy,t)


def ddpm_schedule(steps):
    # Standard linear DDPM beta schedule rescaled to the chosen time discretization.
    beta=torch.linspace(1e-4,.02,steps)*(1000/steps)
    beta=beta.clamp(max=.999);alpha=1-beta
    return beta,alpha,torch.cumprod(alpha,0)


def make_models(kind,data,pack,use_pathway,a):
    first.seed_all(a.seed)
    dims=(pack['gene_dim'],pack['pathway_dim'],pack['conditions'].shape[1],pack['structure'].shape[1],data.target.shape[1])
    model=Generator(dims,use_pathway) if kind=='cgan' else EpsilonDenoiser(dims,use_pathway,a.ddpm_steps)
    opt=torch.optim.AdamW(model.parameters(),lr=.001,weight_decay=.01)
    disc=Discriminator(dims,use_pathway) if kind=='cgan' else None
    dopt=torch.optim.AdamW(disc.parameters(),lr=.001,weight_decay=.01) if disc is not None else None
    return model,opt,disc,dopt


def train_epoch(kind,model,opt,disc,dopt,data,pack,tr,structure,a):
    model.train()
    if disc is not None:disc.train()
    total=0.;_,_,abar=ddpm_schedule(a.ddpm_steps)
    for local in torch.randperm(len(tr)).split(a.batch_size):
        rows=tr[local.numpy()];inputs=first.inputs(data,pack,rows,structure)
        truth=pack['target'][data.target_ids[rows]];weights=pack['weights'][local]
        if kind=='cgan':
            for parameter in disc.parameters():parameter.requires_grad_(True)
            fake=model(*inputs,torch.randn(len(rows),16)).detach()
            real_logits=disc(*inputs,truth);fake_logits=disc(*inputs,fake)
            loss_d=.5*((F.binary_cross_entropy_with_logits(real_logits,torch.ones_like(real_logits),reduction='none')+F.binary_cross_entropy_with_logits(fake_logits,torch.zeros_like(fake_logits),reduction='none'))*weights).mean()
            if not torch.isfinite(loss_d):raise ValueError('Nonfinite discriminator loss')
            dopt.zero_grad();loss_d.backward();nn.utils.clip_grad_norm_(disc.parameters(),1.);dopt.step()
            for parameter in disc.parameters():parameter.requires_grad_(False)
            pred=model(*inputs,torch.randn(len(rows),16));logits=disc(*inputs,pred)
            # Pure conditional GAN: no reconstruction, cycle, module, or pathway loss.
            loss=(F.binary_cross_entropy_with_logits(logits,torch.ones_like(logits),reduction='none')*weights).mean()
        else:
            t=torch.randint(a.ddpm_steps,(len(rows),));noise=torch.randn_like(truth)
            noisy=abar[t,None].sqrt()*truth+(1-abar[t,None]).sqrt()*noise
            epsilon=model(*inputs,noisy,t)
            loss=((epsilon-noise).square().mean(1)*weights).mean()
        if not torch.isfinite(loss):raise ValueError('Nonfinite training loss')
        opt.zero_grad();loss.backward();nn.utils.clip_grad_norm_(model.parameters(),1.);opt.step();total+=float(loss.detach())*len(rows)
    return total/len(tr)


@torch.no_grad()
def draw_mean(kind,model,inputs,a,seed):
    model.eval();gen=torch.Generator().manual_seed(seed)
    context=model.context(*inputs);n=len(context);context=context.repeat_interleave(a.samples,0)
    if kind=='cgan':out=model.decode(context,torch.randn((len(context),16),generator=gen))
    else:
        beta,alpha,abar=ddpm_schedule(a.ddpm_steps)
        out=torch.randn((len(context),model.output_dim),generator=gen)
        for step in reversed(range(a.ddpm_steps)):
            epsilon=model.decode(context,out,torch.full((len(context),),step,dtype=torch.long))
            mean=(out-beta[step]/(1-abar[step]).sqrt()*epsilon)/alpha[step].sqrt()
            out=mean
            if step:
                posterior_variance=beta[step]*(1-abar[step-1])/(1-abar[step])
                out=out+posterior_variance.sqrt()*torch.randn(out.shape,generator=gen)
    if not torch.isfinite(out).all():raise ValueError('Nonfinite DDPM/GAN samples')
    return out.reshape(n,a.samples,-1).mean(1).numpy()


@torch.no_grad()
def predict(kind,model,data,pack,rows,structure,a,seed):
    ids=np.unique(data.target_ids[rows]);lookup={int(cid):i for i,cid in enumerate(ids)}
    sums=np.zeros((len(ids),data.target.shape[1]),float);counts=np.zeros(len(ids),int);pair_sse=0.
    for start in range(0,len(rows),a.batch_size):
        batch=rows[start:start+a.batch_size]
        pred=draw_mean(kind,model,first.inputs(data,pack,batch,structure),a,seed+start)+pack['target_mean']
        targets=data.target_ids[batch];indices=np.array([lookup[int(cid)] for cid in targets])
        np.add.at(sums,indices,pred);np.add.at(counts,indices,1)
        pair_sse+=float(np.sum((pred.astype(float)-data.target[targets])**2))
    return ids,(sums/counts[:,None]).astype(np.float32),float(np.sqrt(pair_sse/(len(rows)*data.target.shape[1])))


def save_json(path,value):
    temp=path.with_suffix('.tmp');temp.write_text(json.dumps(value,indent=2));temp.replace(path)


def load_mlp_reference(a,data,projection):
    root=Path(a.mlp_reference);random_root=Path(a.random_reference)
    cfg=json.loads((root/'sparse'/'config.json').read_text())
    for key in ['prepared','structure','annotation','seed','epochs','check_every','batch_size']:
        if str(cfg[key])!=str(getattr(a,key)):raise ValueError(f'Reference mismatch: {key}')
    if json.loads((root/'pathway_provenance.json').read_text())!=pathway.load_projection(a.annotation,a.prepared)[2]:raise ValueError('Annotation provenance changed')
    random_audit=json.loads((random_root/'randomization_audit.json').read_text())
    if random_audit['original_projection_sha256']!=hashlib.sha256(np.ascontiguousarray(projection).tobytes()).hexdigest():raise ValueError('Random reference annotation mismatch')
    random_mask=np.load(random_root/'random_membership.npy')
    np.testing.assert_array_equal(random_mask.sum(0),(projection>0).sum(0))
    np.testing.assert_array_equal(random_mask.sum(1),(projection>0).sum(1))
    summary=pd.read_csv(random_root/'comparison.csv')
    if set(summary.pathway_mode)!={'none','mean','sparse','random_sparse'}:raise ValueError('Incomplete MLP references')
    reference_split=None
    for mode in ['none','sparse','random_sparse']:
        folder=root/mode if mode!='random_sparse' else random_root/'random'
        split=pd.read_csv(folder/'compound_split.csv').sort_values('compound').reset_index(drop=True)
        if reference_split is None:reference_split=split
        else:pd.testing.assert_frame_equal(reference_split,split)
        filename=f'pathway_{mode}_predictions.npz'
        with np.load(folder/filename) as p:metrics,_=first.evaluate(data,p['target_ids'],p['prediction'])
        record=summary.loc[summary.pathway_mode.eq(mode)].iloc[0]
        for metric in ['absolute_rmse','response_rmse']:np.testing.assert_allclose(metrics[metric],record[metric],rtol=1e-5,atol=1e-6)
    summary=summary.loc[summary.pathway_mode.ne('mean')].copy()
    summary['family']='MLP';summary['origin']='reused_MLP_reference';summary['discriminator_parameters']=0
    return summary,reference_split


def run_matrix(a,families=('cgan','ddpm')):
    global MEMBERSHIP
    if not set(families)<={'cgan','ddpm'}:raise ValueError('Select cgan and/or ddpm')
    projection,pathways,provenance=pathway.load_projection(a.annotation,a.prepared)
    random_mask=np.load(Path(a.random_reference)/'random_membership.npy')
    torch.set_num_threads(a.threads);data=first.Data(a)
    mlp,split=load_mlp_reference(a,data,projection)
    if split.compound.duplicated().any() or set(split.compound)!=set(data.names):raise ValueError('Compound split mismatch')
    def rows(stage):return np.flatnonzero(data.pairs.compound.isin(set(split.loc[split.stage.eq(stage),'compound'])))
    itr,iv,te=rows('inner_train'),rows('inner_validation'),rows('outer_test');tr=np.sort(np.r_[itr,iv])
    if any(set(x)&set(y) for x,y in [(itr,iv),(itr,te),(iv,te)]):raise ValueError('Split overlap')
    a.out.mkdir(parents=True,exist_ok=True)
    config={k:str(v) if isinstance(v,Path) else v for k,v in vars(a).items()}
    config.update(gan='conditional BCE only; no reconstruction loss',diffusion='epsilon-DDPM; scaled linear betas; all reverse steps; no DDIM/v/clipping',denoiser='residual MLP epsilon predictor with identity skip, not post-correction',samples_mean='mean over generated samples per pair, then mean over sources per target')
    config_path=a.out/'config.json'
    if config_path.exists():
        if json.loads(config_path.read_text())!=config:raise ValueError('Resume settings mismatch')
    else:
        if any(a.out.iterdir()):raise ValueError('Unrecognized result directory')
        save_json(config_path,config)
    progress_path=a.out/'progress.json';progress=json.loads(progress_path.read_text()) if progress_path.is_file() else {'runs':{}}
    split.to_csv(a.out/'compound_split.csv',index=False)
    inner=pathway.prepare_pathways(data,itr,iv,projection,pathways,provenance)
    outer=pathway.prepare_pathways(data,tr,te,projection,pathways,provenance)
    joblib.dump(outer['preprocess'],a.out/'preprocess.joblib')
    def collect():
        frames=[mlp]
        for entry in progress['runs'].values():
            path=Path(entry['path'])/'metrics.json'
            if path.is_file():frames.append(pd.DataFrame([json.loads(path.read_text())]))
        table=pd.concat(frames,ignore_index=True);table.to_csv(a.out/'comparison.csv',index=False);return table
    print('Three MLP references verified; no MLP retraining. Structure fixed ON.',flush=True)
    print(f'New families={list(families)}; DDPM full reverse steps={a.ddpm_steps}; generated samples/pair={a.samples}',flush=True)
    collect()
    for kind in families:
        for structure in [True]:
            for mode in ['none','sparse','random_sparse']:
                MEMBERSHIP=random_mask if mode=='random_sparse' else projection
                pathflag='none' if mode=='none' else 'sparse'
                name=f'{kind}_{mode}';entry=progress['runs'].get(name)
                if entry and (Path(entry['path'])/'metrics.json').is_file():print(name,': complete; skipped',flush=True);continue
                if entry:print(name,': partial attempt retained; restarting this variant',flush=True)
                folder=a.out/f'{name}_{uuid.uuid4().hex[:6]}';folder.mkdir()
                progress['runs'][name]={'path':str(folder),'status':'running'};save_json(progress_path,progress)
                print(f'\n{name}: training starts; {folder}',flush=True)
                model,opt,disc,dopt=make_models(kind,data,inner,pathflag,a)
                best=float('inf');epoch_best=None;state=None;history=[]
                for epoch in range(1,a.epochs+1):
                    loss=train_epoch(kind,model,opt,disc,dopt,data,inner,itr,structure,a)
                    if epoch==1 or epoch%a.check_every==0 or epoch==a.epochs:
                        ids,pred,_=predict(kind,model,data,inner,iv,structure,a,a.seed+20000)
                        mse=float(np.mean((pred.astype(float)-data.target[ids])**2))
                        if mse<best:best=mse;epoch_best=epoch;state=copy.deepcopy(model.state_dict())
                        history.append(dict(epoch=epoch,loss=loss,validation_rmse=np.sqrt(mse),best_epoch=epoch_best))
                        pd.DataFrame(history).to_csv(folder/'validation_history.csv',index=False)
                        print(f'{name} epoch={epoch}/{a.epochs} val_RMSE={np.sqrt(mse):.5f} best_epoch={epoch_best}',flush=True)
                torch.save(state,folder/'inner_best.pt')
                model,opt,disc,dopt=make_models(kind,data,outer,pathflag,a)
                for epoch in range(1,epoch_best+1):
                    train_epoch(kind,model,opt,disc,dopt,data,outer,tr,structure,a)
                    if epoch==1 or epoch%10==0 or epoch==epoch_best:print(f'{name} refit={epoch}/{epoch_best}',flush=True)
                ids,pred,pair_rmse=predict(kind,model,data,outer,te,structure,a,a.seed+30000)
                metrics,detail=first.evaluate(data,ids,pred)
                metrics.update(model=name,family='cGAN' if kind=='cgan' else 'epsilon-DDPM',structure=structure,pathway_mode=mode,origin='new_vanilla_run',selected_epoch=epoch_best,validation_rmse=float(np.sqrt(best)),pair_rmse=pair_rmse,parameters=sum(p.numel() for p in model.parameters()),discriminator_parameters=sum(p.numel() for p in disc.parameters()) if disc is not None else 0)
                torch.save(model.state_dict(),folder/'model.pt')
                if disc is not None:torch.save(disc.state_dict(),folder/'discriminator.pt')
                np.savez_compressed(folder/'predictions.npz',target_ids=ids,prediction=pred)
                detail.to_csv(folder/'response_condition_metrics.csv',index=False)
                save_json(folder/'metrics.json',metrics);progress['runs'][name]['status']='complete';save_json(progress_path,progress)
                collect();print(f'{name}: TEST response RMSE={metrics["response_rmse"]:.5f}',flush=True)
    result=collect();print(result.to_string(index=False));print('Saved:',a.out,flush=True);return result


