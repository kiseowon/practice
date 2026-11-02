"""Sparse multi-unit pathway representations, contrasted with mean-only and zero inputs."""
import argparse
import json
from pathlib import Path
from types import SimpleNamespace
import uuid
import numpy as np
import pandas as pd
import torch
from torch import nn
from threadpoolctl import threadpool_limits
import first_structure_mlp as first
import pathway_structure_mlp as pathway


class PatternEncoder(nn.Module):
    def __init__(self,membership,mode,units=4,output_dim=32):
        super().__init__()
        if mode not in ['none','mean','sparse']:raise ValueError(mode)
        edges=np.argwhere(np.asarray(membership)>0)
        if not len(edges) or (np.bincount(edges[:,1],minlength=membership.shape[1])==0).any():raise ValueError('Empty pathway')
        self.mode=mode;self.n_pathways=membership.shape[1];self.units=units
        self.register_buffer('probe_index',torch.tensor(edges[:,0],dtype=torch.long))
        self.register_buffer('pathway_index',torch.tensor(edges[:,1],dtype=torch.long))
        fan_in=np.bincount(edges[:,1],minlength=self.n_pathways)
        weights=torch.randn(len(edges),units)/torch.tensor(np.sqrt(fan_in[edges[:,1]]),dtype=torch.float32)[:,None]
        self.edge_weights=nn.Parameter(weights)
        self.bias=nn.Parameter(torch.zeros(self.n_pathways,units))
        self.readout=nn.Sequential(nn.Dropout(.1),nn.Linear(self.n_pathways*units,output_dim,bias=False),nn.SiLU())

    def pathway_units(self,genes,mean_scores):
        batch=len(genes)
        if self.mode=='none':return genes.new_zeros(batch,self.n_pathways,self.units)
        # Each edge has its own signed, trainable weight. Only KEGG membership edges exist.
        values=genes[:,self.probe_index] if self.mode=='sparse' else mean_scores[:,self.pathway_index]
        weighted=values[:,:,None]*self.edge_weights[None,:,:]
        aggregated=genes.new_zeros(batch,self.n_pathways,self.units).index_add(1,self.pathway_index,weighted)
        return F.silu(aggregated+self.bias)

    def forward(self,genes,mean_scores):return self.readout(self.pathway_units(genes,mean_scores).flatten(1))


# Explicit functional import keeps the standalone notebook and module equivalent.
import torch.nn.functional as F


class PatternMLP(nn.Module):
    def __init__(self,gene_dim,condition_dim,structure_dim,output_dim,membership,mode):
        super().__init__();self.gene_dim=gene_dim
        self.expression_encoder=nn.Sequential(nn.Linear(gene_dim,128),nn.SiLU(),nn.Dropout(.1),nn.Linear(128,64),nn.SiLU())
        self.structure_encoder=nn.Sequential(nn.Linear(structure_dim,64,bias=False),nn.SiLU(),nn.Dropout(.1),nn.Linear(64,32,bias=False),nn.SiLU())
        self.pathway_encoder=PatternEncoder(membership,mode,units=4,output_dim=32)
        self.decoder=nn.Sequential(nn.Linear(64+32+32+condition_dim,64),nn.SiLU(),nn.Dropout(.1),nn.Linear(64,output_dim))
    def forward(self,source,condition,structure):
        genes=source[:,:self.gene_dim];mean_scores=source[:,self.gene_dim:]
        expression=self.expression_encoder(genes);chemical=self.structure_encoder(structure)
        biological=self.pathway_encoder(genes,mean_scores)
        return self.decoder(torch.cat([expression,chemical,biological,condition],1))


def run_patterns(a):
    projection,pathways,provenance=pathway.load_projection(a.annotation,a.prepared)
    if a.out.exists() and any(a.out.iterdir()):raise ValueError('Use a fresh output directory')
    a.out.mkdir(parents=True,exist_ok=True)
    config={k:str(v) if isinstance(v,Path) else v for k,v in vars(a).items()}
    config.update(structure='enabled for all three models',units_per_pathway=4,pathway_readout=32,
        mean_control='same sparse architecture but every member receives its pathway mean score',
        sparse_input='individual TRAIN-standardized probe expression at membership edges',
        scope='PNET-inspired sparse representation, not a reproduction of P-NET/GEARS')
    (a.out/'config.json').write_text(json.dumps(config,indent=2))
    pathways.to_csv(a.out/'pathway_audit.csv',index=False)
    (a.out/'pathway_provenance.json').write_text(json.dumps(provenance,ensure_ascii=False,indent=2))
    print(f'Pathways={len(pathways)}, units/pathway=4, membership edges={int(np.count_nonzero(projection))}',flush=True)
    tables=[];reference_split=None
    for mode in ['none','mean','sparse']:
        folder=a.out/mode
        values={**vars(a),'out':folder,'pathway_mode':mode,'architecture':'matched_sparse_pathway_units4_readout32_structure_on'}
        group=SimpleNamespace(**values)
        def prep(data,tr,ev):return pathway.prepare_pathways(data,tr,ev,projection,pathways,provenance)
        def factory(data,pack,seed):
            first.seed_all(seed)
            return PatternMLP(data.source.shape[1],pack['conditions'].shape[1],pack['structure'].shape[1],data.target.shape[1],projection,mode)
        print(f'\nPathway mode={mode}, structure fixed ON',flush=True)
        result=first.run(group,model_factory=factory,prepare_fn=prep,variants=[(f'pathway_{mode}',True)])
        result['pathway_mode']=mode;result['structure']=True;tables.append(result)
        split=pd.read_csv(folder/'compound_split.csv').sort_values('compound').reset_index(drop=True)
        if reference_split is None:reference_split=split
        else:pd.testing.assert_frame_equal(reference_split,split)
        pd.concat(tables,ignore_index=True).to_csv(a.out/'comparison.csv',index=False)
    result=pd.concat(tables,ignore_index=True)
    if result.parameters.nunique()!=1:raise ValueError('Parameter counts differ')
    print(result.to_string(index=False));print('Saved:',a.out,flush=True)
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for field in ['prepared','structure','annotation','out']:p.add_argument('--'+field,type=Path,required=True)
    p.add_argument('--epochs',type=int,default=30);p.add_argument('--check-every',type=int,default=5);p.add_argument('--batch-size',type=int,default=128);p.add_argument('--threads',type=int,default=4);p.add_argument('--seed',type=int,default=42)
    a=p.parse_args()
    if min(a.epochs,a.check_every,a.batch_size,a.threads)<1:p.error('Positive settings required')
    with threadpool_limits(limits=a.threads):run_patterns(a)
