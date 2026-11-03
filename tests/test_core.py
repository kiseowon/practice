import json
from pathlib import Path
import sys
import unittest
import numpy as np
import torch
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from pathway_pattern_mlp import PatternEncoder
from pathway_random_control import shuffle_membership

class CoreTests(unittest.TestCase):
    def test_gene_patterns_and_membership(self):
        m=np.array([[1.],[1.],[0.]])
        encoder=PatternEncoder(m,'sparse',units=2)
        with torch.no_grad():encoder.edge_weights.copy_(torch.eye(2))
        x=torch.tensor([[2.,-2.,0.],[-2.,2.,10.]])
        scores=torch.zeros(2,1)
        out=encoder.pathway_units(x,scores)
        self.assertFalse(torch.allclose(out[0],out[1]))
        changed=x.clone();changed[:,2]+=100
        self.assertTrue(torch.equal(out,encoder.pathway_units(changed,scores)))
        out.sum().backward()
        self.assertTrue(torch.isfinite(encoder.edge_weights.grad).all())
    def test_randomization_degrees(self):
        original=np.random.default_rng(17).random((30,12))<.2
        randomized,_=shuffle_membership(original,seed=4,swaps_per_edge=2)
        np.testing.assert_array_equal(original.sum(0),randomized.sum(0))
        np.testing.assert_array_equal(original.sum(1),randomized.sum(1))
        self.assertFalse(np.array_equal(original,randomized))
    def test_standalone_notebook_definitions(self):
        # Archived notebooks contain historical training/API calls; inspect only
        # the three standalone definition notebooks this test was written for.
        for name in ['pathway_pattern_mlp.ipynb', 'pathway_random_control.ipynb',
                     'sparse_pathway_gan_ddpm.ipynb']:
            file=ROOT/'notebooks'/name
            nb=json.loads(file.read_text());namespace={'__name__':'test_notebook'}
            code=[c for c in nb['cells'] if c['cell_type']=='code']
            for c in code[:-1]:exec(compile(''.join(c['source']),str(file),'exec'),namespace)
            compile(''.join(code[-1]['source']),str(file),'exec')

if __name__=='__main__':unittest.main()
