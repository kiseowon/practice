"""기존 잠재 diffusion과 같은 데이터·설정 선택으로 작은 MLP 또는 혼합손실 cGAN을 비교합니다."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import torch
from torch import nn
import torch.nn.functional as F
import train_conditioned_diffusion as experiment

MODEL_KIND = 'mlp'


class SmallModel(nn.Module):
    def __init__(self, condition_dim, latent_dim, steps):
        super().__init__()
        self.noise_dim = 16 if MODEL_KIND == 'cgan' else 0
        self.generator = nn.Sequential(
            nn.Linear(condition_dim+self.noise_dim,128),nn.SiLU(),nn.Dropout(.1),
            nn.Linear(128,64),nn.SiLU(),nn.Linear(64,latent_dim))
        if self.noise_dim:
            # 판별기는 생성기 optimizer와 분리해 관리합니다.
            discriminator = nn.Sequential(nn.Linear(condition_dim+latent_dim,64),nn.LeakyReLU(.2),nn.Linear(64,32),nn.LeakyReLU(.2),nn.Linear(32,1))
            object.__setattr__(self,'discriminator',discriminator)
            object.__setattr__(self,'discriminator_optimizer',None)

    def generate(self, c, noise=None):
        if self.noise_dim:
            if noise is None:
                noise=torch.randn((len(c),self.noise_dim),device=c.device)
            c=torch.cat([c,noise],1)
        return self.generator(c)


def fit_epoch(model,optimizer,c,y,weights,steps):
    model.train()
    if model.noise_dim:
        d=model.discriminator.to(c.device)
        d.train()
        if model.discriminator_optimizer is None:
            object.__setattr__(model,'discriminator_optimizer',torch.optim.AdamW(d.parameters(),lr=.001,weight_decay=.01))
    total=0.
    for idx in torch.randperm(len(c),device=c.device).split(64):
        cb,yb,wb=c[idx],y[idx],weights[idx]
        if model.noise_dim:
            d=model.discriminator
            for p in d.parameters():
                p.requires_grad_(True)
            fake=model.generate(cb).detach()
            real_logits=d(torch.cat([cb,yb],1)).squeeze(1)
            fake_logits=d(torch.cat([cb,fake],1)).squeeze(1)
            dloss=(wb*(F.binary_cross_entropy_with_logits(real_logits,torch.ones_like(real_logits),reduction='none')+F.binary_cross_entropy_with_logits(fake_logits,torch.zeros_like(fake_logits),reduction='none'))*.5).mean()
            model.discriminator_optimizer.zero_grad()
            dloss.backward()
            nn.utils.clip_grad_norm_(d.parameters(),1.)
            model.discriminator_optimizer.step()
            for p in d.parameters():
                p.requires_grad_(False)
        prediction=model.generate(cb)
        reconstruction=(prediction-yb).square().mean(1)
        loss=(wb*reconstruction).mean()
        if model.noise_dim:
            logits=d(torch.cat([cb,prediction],1)).squeeze(1)
            adversarial=F.binary_cross_entropy_with_logits(logits,torch.ones_like(logits),reduction='none')
            # 점예측 비교를 위한 재구성 중심의 고정 혼합손실
            loss=loss+.1*(wb*adversarial).mean()
        optimizer.zero_grad()
        loss.backward()
        nn.utils.clip_grad_norm_(model.parameters(),1.)
        optimizer.step()
        total+=loss.item()*len(idx)
    return total/len(c)


@torch.no_grad()
def sample(model,condition,latent_dim,steps,samples,seed,device):
    model.eval()
    generator=torch.Generator(device=device).manual_seed(seed)
    output=[]
    for start in range(0,len(condition),128):
        c=torch.tensor(condition[start:start+128],dtype=torch.float32,device=device)
        n=len(c)
        if model.noise_dim:
            c=c.repeat_interleave(samples,0)
            noise=torch.randn((len(c),model.noise_dim),generator=generator,device=device)
            pred=model.generate(c,noise).reshape(n,samples,latent_dim).mean(1)
        else:
            pred=model.generate(c)
        if not torch.isfinite(pred).all():
            raise ValueError('모델 예측 NaN/inf')
        output.append(pred.cpu().numpy())
    return np.concatenate(output)


def main():
    global MODEL_KIND
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--prepared',type=Path,required=True)
    p.add_argument('--structure',type=Path,required=True)
    p.add_argument('--structure-kind',choices=['binary','count'],default='binary')
    p.add_argument('--no-structure',action='store_true')
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--model',choices=['mlp','cgan'],required=True)
    p.add_argument('--seed',type=int,default=42)
    p.add_argument('--model-seed',type=int,default=42)
    p.add_argument('--folds',type=int,default=5)
    p.add_argument('--input-components',type=int,default=64)
    p.add_argument('--latent-components',type=int,default=32)
    p.add_argument('--epochs',type=int,default=100)
    p.add_argument('--check-every',type=int,default=25)
    p.add_argument('--samples',type=int,default=8)
    p.add_argument('--threads',type=int,default=4)
    p.add_argument('--device',choices=['cpu','cuda'],default='cpu')
    p.add_argument('--alphas',type=float,nargs='+',default=[1.,10.,100.,1000.])
    p.add_argument('--strengths',type=float,nargs='+',default=[0.,.1,.25,.5,1.])
    a=p.parse_args();a.steps=30
    if min(a.folds,a.steps)<2 or min(a.input_components,a.latent_components,a.epochs,a.check_every,a.samples,a.threads)<1 or min(a.alphas)<=0 or min(a.strengths)<0:
        p.error('설정 범위 오류')
    if a.out.exists() and any(a.out.iterdir()):
        p.error('기존 결과 보호를 위해 비어 있는 출력 폴더를 사용하세요.')
    MODEL_KIND=a.model
    a.neural_name=f'latent_{MODEL_KIND}'
    experiment.Denoiser=SmallModel
    experiment.fit_epoch=fit_epoch
    experiment.sample=sample
    print('선택 모델:',MODEL_KIND,flush=True)
    experiment.run(a)
    name=f'latent_{MODEL_KIND}'
    for path in a.out.glob('fold_*_latent_diffusion.npz'):
        path.rename(path.with_name(path.name.replace('latent_diffusion',name)))
    metrics=pd.read_csv(a.out/'fold_metrics.csv')
    metrics['model']=metrics.model.replace({'latent_diffusion':name})
    metrics.to_csv(a.out/'fold_metrics.csv',index=False)
    settings=pd.read_csv(a.out/'selected_settings.csv').rename(columns={'diffusion_epoch':f'{MODEL_KIND}_epoch','diffusion_strength':f'{MODEL_KIND}_strength'})
    settings.to_csv(a.out/'selected_settings.csv',index=False)
    print('저장 모델 이름:',name,flush=True)
    print('cGAN은 재구성 MSE + 0.1 × 조건부 적대 손실이며 원본 AIVIVE GAN의 재현이 아닙니다.')


if __name__=='__main__':
    main()
