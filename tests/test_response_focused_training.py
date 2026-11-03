import contextlib
import io
from pathlib import Path
import sys
import tempfile
import unittest
from types import SimpleNamespace
from unittest.mock import Mock, patch

import numpy as np
import pandas as pd
import torch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import response_focused_training as response
import first_structure_mlp as first


class ResponseProtocolTests(unittest.TestCase):
    def data(self):
        return SimpleNamespace(
            source=np.array([[1., 2.], [3., 4.], [9., 10.]], dtype=np.float32),
            target=np.array([[1., 3.], [5., 7.], [900., 1000.]], dtype=np.float32),
            source_ids=np.array([0, 1, 1, 2]), target_ids=np.array([0, 1, 1, 2]),
            structure=np.zeros((3, 1), dtype=np.float32),
            features=pd.DataFrame(dict(source_dose=[0., 1., 1., 10.], source_time=[2., 8., 8., 24.],
                target_dose=[0., 1., 1., 10.], target_time=[3., 6., 6., 24.],
                source_unit=['u', 'u', 'u', 'new'], target_unit=['v']*4,
                source_level=['Control', 'Low', 'Low', 'High'], target_level=['Control', 'Low', 'Low', 'High'])))

    def test_target_scale_fits_unique_training_targets_only(self):
        d=self.data();projection=np.array([[.5],[.5]]);paths=pd.DataFrame({'pathway':['p']})
        a=response.prepare_stage(d,np.array([0,1,2]),np.array([3]),projection,paths,{})
        self.assertEqual(a['target_scale'],2.)
        np.testing.assert_allclose(a['target'][:2], [[-1.,-1.],[1.,1.]])
        # One common scale preserves relative probe errors; held-out changes cannot fit it.
        d.target[2]=[-999., 999.];d.source[2]=[2000.,3000.]
        b=response.prepare_stage(d,np.array([0,1,2]),np.array([3]),projection,paths,{})
        self.assertEqual(a['target_scale'],b['target_scale'])
        np.testing.assert_array_equal(a['source'][:2],b['source'][:2])
        d.target[:2]=4
        with self.assertRaisesRegex(ValueError,'degenerate'):
            response.prepare_stage(d,np.array([0,1,2]),np.array([3]),projection,paths,{})

    def test_prediction_inverse_scale_and_source_average(self):
        d=SimpleNamespace(target=np.zeros((1,2)),target_ids=np.array([0,0]))
        model=Mock();model.return_value=torch.tensor([[1.,2.],[3.,4.]])
        pack={'target_scale':2.,'target_mean':np.array([10.,20.])}
        with patch.object(first,'inputs',return_value=(None,None,None)):
            ids,pred=response.predict('mlp',model,d,pack,np.array([0,1]),77)
        np.testing.assert_array_equal(ids,[0]);np.testing.assert_allclose(pred,[[14.,26.]])

    def test_selection_uses_response_not_absolute_and_never_test(self):
        model=Mock();model.state_dict.return_value={}
        protocol={**response.DEFAULT_PROTOCOL,'learning_rates':[.001,.0003], 'max_epochs':1}
        validation=np.array([1]);data=SimpleNamespace()
        values=[({'absolute_rmse':.1,'response_rmse':.8},pd.DataFrame()),
                ({'absolute_rmse':.9,'response_rmse':.2},pd.DataFrame())]
        with tempfile.TemporaryDirectory() as directory:
            with patch.dict(response.DEFAULT_PROTOCOL,protocol), \
                 patch.object(response,'make_model',return_value=(model,Mock(),None,None)), \
                 patch.object(response,'train_epoch',return_value=0.) as train, \
                 patch.object(response,'predict',return_value=(np.array([0]),np.zeros((1,2)))) as predict, \
                 patch.object(first,'evaluate',side_effect=values):
                choice=response.select_candidate('mlp',data,{},np.array([0]),validation,
                    np.ones((2,1)),'sparse',42,Path(directory))
        self.assertEqual(choice['learning_rate'],.0003)
        self.assertEqual(train.call_count,2)
        for call in predict.call_args_list:np.testing.assert_array_equal(call.args[4],validation)

    def test_final_refit_restores_selection_and_evaluates_once(self):
        model=Mock();model.state_dict.return_value={}
        choice={'learning_rate':.0003,'selected_epoch':2,'validation_response_rmse':.2,'candidate':0}
        pack={'target_scale':2.}
        with tempfile.TemporaryDirectory() as directory:
            with patch.object(response,'make_model',return_value=(model,Mock(),None,None)) as make, \
                 patch.object(response,'train_epoch',return_value=0.) as train, \
                 patch.object(response,'predict',return_value=(np.array([2]),np.ones((1,2)))) as predict, \
                 patch.object(first,'evaluate',return_value=({'response_rmse':.5,'absolute_rmse':.4},pd.DataFrame())):
                response.final_refit('mlp',None,pack,np.array([0,1]),np.array([2]),
                                    np.ones((2,1)),'sparse',42,choice,Path(directory))
            self.assertEqual(train.call_count,2);predict.assert_called_once()
            self.assertEqual(make.call_args.args[-1],.0003)
            for call in train.call_args_list:np.testing.assert_array_equal(call.args[-1],[0,1])

    def test_default_plan_never_reads_data_or_trains(self):
        with tempfile.TemporaryDirectory() as directory:
            out=Path(directory)/'new'
            with patch.object(response,'execute') as execute, contextlib.redirect_stdout(io.StringIO()):
                planned=response.main(['--prepared','/missing','--structure','/missing.csv',
                    '--annotation','/missing-a','--out',str(out)])
            self.assertEqual(len(planned),15);execute.assert_not_called();self.assertFalse(out.exists())


if __name__=='__main__':unittest.main()
