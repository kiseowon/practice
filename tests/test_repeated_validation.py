import contextlib
import io
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import repeated_pathway_validation as repeat


class RepeatedValidationTests(unittest.TestCase):
    def test_plan_never_reads_inputs_creates_outputs_or_trains(self):
        with tempfile.TemporaryDirectory() as directory:
            out = Path(directory) / 'not-created'
            with patch.object(repeat, 'execute', side_effect=AssertionError('training called')):
                with contextlib.redirect_stdout(io.StringIO()):
                    planned = repeat.main(['--prepared', '/missing/prepared', '--structure', '/missing/s.csv',
                                           '--annotation', '/missing/kegg', '--out', str(out)])
            self.assertEqual(len(planned), 75)
            self.assertFalse(out.exists())
            self.assertEqual(len({tuple(j.items()) for j in planned}), 75)

    def test_baselines_not_duplicated_by_mask(self):
        planned = repeat.jobs([42], [10, 11], [1, 2, 3])
        self.assertEqual(len(planned), 30)
        for family in repeat.FAMILIES:
            for model_seed in [10, 11]:
                subset = [j for j in planned if j['family'] == family and j['model_seed'] == model_seed]
                self.assertEqual([j['pathway_mode'] for j in subset].count('none'), 1)
                self.assertEqual([j['pathway_mode'] for j in subset].count('sparse'), 1)
                self.assertEqual([j['pathway_mode'] for j in subset].count('random_sparse'), 3)

    def test_compound_split_matches_original_and_is_seed_specific(self):
        names = [f'compound-{i:03d}' for i in range(140)]
        split = repeat.compound_split(names[::-1], 42)
        self.assertEqual(split.stage.value_counts().to_dict(),
                         {'inner_train': 89, 'outer_test': 28, 'inner_validation': 23})
        self.assertFalse(split.compound.duplicated().any())
        perm = np.random.default_rng(42).permutation(sorted(names))
        self.assertEqual(set(split.loc[split.stage.eq('outer_test'), 'compound']), set(perm[:28]))
        pd.testing.assert_frame_equal(split, repeat.compound_split(names, 42))
        self.assertNotEqual(set(perm[:28]), set(repeat.compound_split(names, 43)
                                               .query("stage == 'outer_test'").compound))

    def test_summary_pairs_masks_within_split_not_as_independent_baselines(self):
        records = []
        for split, offset in [(42, 0), (43, 10)]:
            for mode, values in [('none', [5]), ('sparse', [3]), ('random_sparse', [2, 4, 9])]:
                for i, value in enumerate(values):
                    records.append(dict(split_seed=split, model_seed=42, family='MLP',
                                        pathway_mode=mode, graph_seed=i if mode == 'random_sparse' else None,
                                        absolute_rmse=value + offset, response_rmse=value + offset))
        unit, paired, split_paired, summary = repeat.summarize(pd.DataFrame(records))
        self.assertEqual(len(unit), 6)
        self.assertEqual(len(paired), 4)
        self.assertEqual(len(split_paired), 4)
        self.assertTrue((summary.absolute_rmse_delta_count == 2).all())
        self.assertAlmostEqual(summary.set_index('reference').loc['none', 'response_rmse_delta_mean'], -2)
        self.assertAlmostEqual(summary.set_index('reference').loc['random_sparse', 'response_rmse_delta_mean'], -2)
        self.assertTrue((summary.response_rmse_delta_std == 0).all())
        with self.assertRaises(ValueError):
            repeat.summarize(pd.DataFrame(records).query("pathway_mode != 'none'"))

    def test_stage_selection_and_refit_use_correct_data_without_training(self):
        from types import SimpleNamespace
        from unittest.mock import Mock
        import first_structure_mlp as first
        import sparse_pathway_information_matrix as gen
        job = dict(split_seed=43, model_seed=77, family='epsilon-DDPM',
                   pathway_mode='random_sparse', graph_seed=4)
        inner, outer = {'stage': 'inner'}, {'stage': 'outer'}
        itr, iv, tr, te = (np.array(x) for x in [[0], [1], [0, 1], [2]])
        model = Mock()
        model.state_dict.return_value = {}
        settings = SimpleNamespace(epochs=2, check_every=1, batch_size=2, samples=4, ddpm_steps=100)
        calls = []
        def predict(kind, network, data, pack, rows, structure, args, seed):
            calls.append((pack['stage'], rows.tolist(), args.seed, seed))
            return np.array([0]), np.zeros((1, 2)), 0.0
        with tempfile.TemporaryDirectory() as directory:
            with patch.object(gen, 'make_models', return_value=(model, Mock(), None, None)) as make, \
                 patch.object(gen, 'train_epoch', return_value=0.0) as epoch, \
                 patch.object(gen, 'predict', side_effect=predict), \
                 patch.object(first, 'evaluate', return_value=(dict(absolute_rmse=0., response_rmse=0.), pd.DataFrame())):
                record = repeat.fit_job(job, SimpleNamespace(target=np.zeros((1, 2))), inner, outer,
                                        (itr, iv, tr, te), np.ones((2, 1)), settings, Path(directory))
            self.assertEqual(record['selected_epoch'], 1)
            self.assertEqual(calls, [('inner', [1], 77, 20077), ('inner', [1], 77, 20077),
                                     ('outer', [2], 77, 30077)])
            self.assertEqual(make.call_count, 2)
            self.assertEqual([c.args[7].tolist() for c in epoch.call_args_list], [[0], [0], [0, 1]])
            self.assertEqual(epoch.call_count, 3)
            self.assertTrue((Path(directory) / 'metrics.json').is_file())

    def test_output_overlap_and_foreign_output_rejected_before_training(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            prepared, annotation, output = (root / n for n in ['prepared', 'annotation', 'out'])
            prepared.mkdir(); annotation.mkdir(); output.mkdir()
            structure = root / 'structure.csv'; structure.write_text('fixture')
            for name in ['vitro_expression.npy', 'vivo_expression.npy', 'vitro_conditions.csv',
                         'vivo_conditions.csv', 'pairs.csv.gz', 'probes.csv']:
                (prepared / name).write_text('fixture')
            for name in ['projection.npy', 'pathway_audit.csv', 'provenance.json']:
                (annotation / name).write_text('fixture')
            args = repeat.parser().parse_args(['--prepared', str(prepared), '--structure', str(structure),
                                               '--annotation', str(annotation), '--out', str(prepared / 'child')])
            with self.assertRaisesRegex(ValueError, 'separate'):
                repeat.execute(args, [])
            args.out = output
            (output / 'existing.txt').write_text('preserve me')
            with self.assertRaisesRegex(ValueError, 'different configuration'):
                repeat.execute(args, [])
            self.assertEqual((output / 'existing.txt').read_text(), 'preserve me')
            self.assertFalse((output / 'run.json').exists())

    def test_duplicate_seeds_rejected_before_execute(self):
        with patch.object(repeat, 'execute', side_effect=AssertionError('called')):
            with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
                repeat.main(['--prepared', '/p', '--structure', '/s', '--annotation', '/a',
                             '--out', '/o', '--split-seeds', '42', '42', '--execute'])


if __name__ == '__main__':
    unittest.main()
