"""Degree-preserving random pathway control for the matched sparse MLP."""
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import pandas as pd
import first_structure_mlp as first
import pathway_structure_mlp as pathway
from pathway_pattern_mlp import PatternMLP


def shuffle_membership(projection,seed=20261102,swaps_per_edge=10):
    original=np.asarray(projection)>0
    mask=original.copy();edges=np.argwhere(mask)
    rng=np.random.default_rng(seed);accepted=0;attempts=0
    goal=int(len(edges)*swaps_per_edge);limit=max(1000,goal*100)
    while accepted<goal and attempts<limit:
        attempts+=1
        i,j=rng.integers(len(edges),size=2)
        a,b=edges[i];c,d=edges[j]
        if a==c or b==d or mask[a,d] or mask[c,b]:continue
        mask[a,b]=mask[c,d]=False;mask[a,d]=mask[c,b]=True
        edges[i]=[a,d];edges[j]=[c,b];accepted+=1
        if accepted%max(1,goal//10)==0:print(f'Randomization: {accepted}/{goal} accepted swaps',flush=True)
    if accepted<goal:raise RuntimeError('Not enough accepted swaps; no training started')
    np.testing.assert_array_equal(mask.sum(0),original.sum(0))
    np.testing.assert_array_equal(mask.sum(1),original.sum(1))
    overlap=int((mask&original).sum())
    if overlap==int(original.sum()):raise RuntimeError('Randomization did not change membership')
    audit=dict(seed=seed,accepted_swaps=accepted,attempts=attempts,n_edges=len(edges),
        original_edge_overlap=overlap/len(edges),probe_degrees_preserved=True,pathway_degrees_preserved=True,
        note='Degree-preserving double-edge swaps; overlap can remain; not guaranteed uniform sampling')
    return mask.astype(np.float64),audit


def run_random_control(a):
    previous=Path(a.previous)
    old=json.loads((previous/'sparse'/'config.json').read_text())
    for field in ['prepared','structure','annotation','seed','epochs','check_every','batch_size']:
        if str(old[field])!=str(getattr(a,field)):raise ValueError(f'Reference configuration differs: {field}')
    baseline=pd.read_csv(previous/'comparison.csv')
    if set(baseline.pathway_mode)!={'none','mean','sparse'} or len(baseline)!=3:raise ValueError('Expected three original reference models')
    projection,pathways,provenance=pathway.load_projection(a.annotation,a.prepared)
    saved_provenance=json.loads((previous/'pathway_provenance.json').read_text())
    if provenance!=saved_provenance:raise ValueError('Reference annotation provenance differs')
    if a.out.exists() and any(a.out.iterdir()):raise ValueError('Use a new output folder')
    a.out.mkdir(parents=True,exist_ok=True)
    randomized,audit=shuffle_membership(projection,a.graph_seed)
    audit['original_projection_sha256']=hashlib.sha256(np.ascontiguousarray(projection).tobytes()).hexdigest()
    np.save(a.out/'random_membership.npy',randomized.astype(np.uint8))
    (a.out/'randomization_audit.json').write_text(json.dumps(audit,indent=2))
    (a.out/'config.json').write_text(json.dumps({k:str(v) if isinstance(v,Path) else v for k,v in vars(a).items()},indent=2))
    group=SimpleNamespace(**{**vars(a),'out':a.out/'random','pathway_mode':'random_sparse'})
    # Original pathway means are prepared for compatibility, but sparse mode never reads them.
    def prep(data,tr,ev):return pathway.prepare_pathways(data,tr,ev,projection,pathways,provenance)
    def factory(data,pack,seed):
        first.seed_all(seed)
        return PatternMLP(data.source.shape[1],pack['conditions'].shape[1],pack['structure'].shape[1],data.target.shape[1],randomized,'sparse')
    # Check the reference split BEFORE training as well as checking emitted files afterwards.
    data=first.Data(group)
    names=np.random.default_rng(a.seed).permutation(data.names)
    n_test=max(1,int(np.ceil(.2*len(names))))
    remaining=names[n_test:];n_val=max(1,int(np.ceil(.2*len(remaining))))
    sets=[set(remaining[n_val:]),set(remaining[:n_val]),set(names[:n_test])]
    expected=pd.DataFrame([dict(compound=n,stage=stage) for stage,values in zip(['inner_train','inner_validation','outer_test'],sets) for n in sorted(values)])
    reference=pd.read_csv(previous/'sparse'/'compound_split.csv').sort_values('compound').reset_index(drop=True)
    pd.testing.assert_frame_equal(reference,expected.sort_values('compound').reset_index(drop=True))
    print('Training random sparse model only; original three models are reused.',flush=True)
    result=first.run(group,model_factory=factory,prepare_fn=prep,variants=[('pathway_random_sparse',True)])
    reference=pd.read_csv(previous/'sparse'/'compound_split.csv').sort_values('compound').reset_index(drop=True)
    emitted=pd.read_csv(group.out/'compound_split.csv').sort_values('compound').reset_index(drop=True)
    pd.testing.assert_frame_equal(reference,emitted)
    result['pathway_mode']='random_sparse';result['structure']=True
    combined=pd.concat([baseline,result],ignore_index=True)
    if combined.parameters.nunique()!=1:raise ValueError('Parameter counts differ')
    combined.to_csv(a.out/'comparison.csv',index=False)
    print(combined[['pathway_mode','absolute_rmse','response_rmse','selected_epoch']].to_string(index=False),flush=True)
    print('Edge overlap with KEGG:',audit['original_edge_overlap'],flush=True)
    print('Saved:',a.out,flush=True)
    return combined
