"""New response-focused protocol; default plans only, never runs research training.

Reuses established encoders/losses while making target scale and validation
selection explicit. This is a NEW experiment, not reproduction of old results.
"""
import argparse
import json
from pathlib import Path
from types import SimpleNamespace
import uuid

LABELS = {'mlp': 'MLP', 'cgan': 'cGAN', 'ddpm': 'epsilon-DDPM'}
DEFAULT_PROTOCOL = {
    'learning_rates': [0.0003, 0.001],
    'dropout': 0.1, 'weight_decay': 0.01, 'batch_size': 128,
    'max_epochs': 30, 'check_every': 5, 'gradient_clip': 1,
    'pathway_units': 4, 'pathway_readout': 32,
    'samples': 4, 'ddpm_steps': 100,
    'target_transform': 'training_mean_then_shared_rms',
    'selection': 'validation_response_rmse',
}


def prepare_stage(data, training, evaluation, projection, pathways, provenance):
    import numpy as np
    import pathway_structure_mlp as pathway
    pack = pathway.prepare_pathways(data, training, evaluation, projection, pathways, provenance)
    ids = np.unique(data.target_ids[training])
    residual = data.target[ids].astype(float) - pack['target_mean'].astype(float)
    scale = float(np.sqrt(np.mean(residual ** 2)))
    if not np.isfinite(scale) or scale <= 1e-8:
        raise ValueError('Training target RMS is degenerate; inspect data before proceeding')
    pack['target'] = pack['target'] / scale
    pack['target_scale'] = scale
    pack['preprocess']['target_scale'] = scale
    pack['preprocess']['target_transform'] = DEFAULT_PROTOCOL['target_transform']
    pack['preprocess']['source_constant_probe_count'] = int(
        np.count_nonzero(pack['preprocess']['source_scale'].var_ == 0))
    return pack


def make_model(kind, data, pack, membership, mode, seed, lr):
    import torch
    import first_structure_mlp as first
    import sparse_pathway_information_matrix as generative
    from pathway_pattern_mlp import PatternMLP
    settings = SimpleNamespace(seed=seed, ddpm_steps=DEFAULT_PROTOCOL['ddpm_steps'])
    if kind == 'mlp':
        first.seed_all(seed)
        model = PatternMLP(pack['gene_dim'], pack['conditions'].shape[1],
                           pack['structure'].shape[1], data.target.shape[1], membership, mode)
        discriminator = None
    else:
        generative.MEMBERSHIP = membership
        model, _, discriminator, _ = generative.make_models(kind, data, pack, mode, settings)
    for network in [model, discriminator]:
        if network is not None:
            for layer in network.modules():
                if isinstance(layer, torch.nn.Dropout):
                    layer.p = DEFAULT_PROTOCOL['dropout']
    opt = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=DEFAULT_PROTOCOL['weight_decay'])
    dopt = torch.optim.AdamW(discriminator.parameters(), lr=lr,
                            weight_decay=DEFAULT_PROTOCOL['weight_decay']) if discriminator is not None else None
    return model, opt, discriminator, dopt


def train_epoch(kind, model, opt, disc, dopt, data, pack, rows):
    import first_structure_mlp as first
    import sparse_pathway_information_matrix as generative
    if kind == 'mlp':
        return first.train_epoch(model, opt, data, pack, rows, True, DEFAULT_PROTOCOL['batch_size'])
    settings = SimpleNamespace(batch_size=DEFAULT_PROTOCOL['batch_size'], ddpm_steps=DEFAULT_PROTOCOL['ddpm_steps'])
    return generative.train_epoch(kind, model, opt, disc, dopt, data, pack, rows, True, settings)


def predict(kind, model, data, pack, rows, sample_seed):
    import numpy as np
    import torch
    import first_structure_mlp as first
    import sparse_pathway_information_matrix as generative
    ids = np.unique(data.target_ids[rows])
    lookup = {int(cid): i for i, cid in enumerate(ids)}
    sums = np.zeros((len(ids), data.target.shape[1]), dtype=float)
    counts = np.zeros(len(ids), dtype=int)
    settings = SimpleNamespace(samples=DEFAULT_PROTOCOL['samples'], ddpm_steps=DEFAULT_PROTOCOL['ddpm_steps'])
    model.eval()
    with torch.no_grad():
        for start in range(0, len(rows), DEFAULT_PROTOCOL['batch_size']):
            batch = rows[start:start + DEFAULT_PROTOCOL['batch_size']]
            inputs = first.inputs(data, pack, batch, True)
            normalized = model(*inputs).numpy() if kind == 'mlp' else generative.draw_mean(
                kind, model, inputs, settings, sample_seed + start)
            restored = normalized * pack['target_scale'] + pack['target_mean']
            if not np.isfinite(restored).all():
                raise ValueError('Nonfinite restored predictions')
            positions = np.array([lookup[int(cid)] for cid in data.target_ids[batch]])
            np.add.at(sums, positions, restored)
            np.add.at(counts, positions, 1)
    return ids, (sums / counts[:, None]).astype(np.float32)


def select_candidate(kind, data, pack, training, validation, membership, mode, seed, folder):
    """No outer test IDs are passed here. Joint selection of LR and epoch on inner validation."""
    import pandas as pd
    import torch
    import first_structure_mlp as first
    import numpy as np
    history, best = [], None
    for index, lr in enumerate(DEFAULT_PROTOCOL['learning_rates']):
        model, opt, disc, dopt = make_model(kind, data, pack, membership, mode, seed, lr)
        for epoch in range(1, DEFAULT_PROTOCOL['max_epochs'] + 1):
            loss = train_epoch(kind, model, opt, disc, dopt, data, pack, training)
            if epoch == 1 or epoch % DEFAULT_PROTOCOL['check_every'] == 0 or epoch == DEFAULT_PROTOCOL['max_epochs']:
                ids, prediction = predict(kind, model, data, pack, validation, seed + 20000)
                metrics, _ = first.evaluate(data, ids, prediction)
                score = metrics['response_rmse']
                if not np.isfinite(score):
                    raise ValueError('Nonfinite validation response RMSE')
                history.append(dict(candidate=index, learning_rate=lr, epoch=epoch, loss=loss, **metrics))
                if best is None or score < best['validation_response_rmse']:
                    best = dict(learning_rate=lr, selected_epoch=epoch,
                                validation_response_rmse=score, candidate=index)
                    torch.save(model.state_dict(), folder / 'inner_best.pt')
                pd.DataFrame(history).to_csv(folder / 'validation_history.csv', index=False)
    return best


def final_refit(kind, data, pack, training, testing, membership, mode, seed, choice, folder):
    import numpy as np
    import torch
    import first_structure_mlp as first
    model, opt, disc, dopt = make_model(kind, data, pack, membership, mode, seed, choice['learning_rate'])
    for epoch in range(choice['selected_epoch']):
        train_epoch(kind, model, opt, disc, dopt, data, pack, training)
    ids, prediction = predict(kind, model, data, pack, testing, seed + 30000)
    metrics, detail = first.evaluate(data, ids, prediction)
    np.savez_compressed(folder / 'predictions.npz', target_ids=ids, prediction=prediction)
    detail.to_csv(folder / 'response_condition_metrics.csv', index=False)
    torch.save(model.state_dict(), folder / 'model.pt')
    if disc is not None:
        torch.save(disc.state_dict(), folder / 'discriminator.pt')
    return dict(**metrics, **choice, target_scale=pack['target_scale'])


def execute(a, planned):
    import numpy as np
    import pandas as pd
    import torch
    import sklearn
    import joblib
    import platform
    import first_structure_mlp as first
    import pathway_structure_mlp as pathway
    from pathway_random_control import shuffle_membership
    import repeated_pathway_validation as repeated
    from threadpoolctl import threadpool_limits
    out = a.out.resolve()
    for source in [a.prepared, a.structure, a.annotation]:
        source = source.resolve()
        if source == out or source in out.parents or out in source.parents:
            raise ValueError('Output must be separate from inputs')
    inputs = [a.prepared / name for name in ['vitro_expression.npy', 'vivo_expression.npy',
              'vitro_conditions.csv', 'vivo_conditions.csv', 'pairs.csv.gz', 'probes.csv']]
    inputs += [a.structure] + [a.annotation / name for name in
                              ['projection.npy', 'pathway_audit.csv', 'provenance.json']]
    modules = ['response_focused_training.py', 'repeated_pathway_validation.py', 'first_structure_mlp.py',
               'pathway_structure_mlp.py', 'pathway_pattern_mlp.py', 'sparse_pathway_information_matrix.py',
               'pathway_random_control.py']
    files = inputs + [Path(__file__).with_name(name) for name in modules]
    contract = dict(protocol=DEFAULT_PROTOCOL,
                    arguments={k: str(v.resolve()) if isinstance(v, Path) else v
                               for k, v in vars(a).items() if k != 'execute'},
                    sha256={str(p.resolve()): repeated.digest(p) for p in files},
                    versions=dict(python=platform.python_version(), numpy=np.__version__, pandas=pd.__version__,
                                  torch=torch.__version__, sklearn=sklearn.__version__, joblib=joblib.__version__))
    if out.exists() and any(out.iterdir()):
        if not (out / 'run.json').is_file() or json.loads((out / 'run.json').read_text()) != contract:
            raise ValueError('Output contract differs; use a new directory')
    else:
        out.mkdir(parents=True, exist_ok=True)
        repeated.save_json(out / 'run.json', contract)
    data = first.Data(a)
    projection, pathways, provenance = pathway.load_projection(a.annotation, a.prepared)
    masks = {}
    for graph in a.graph_seeds:
        file = out / f'mask_{graph}.npy'
        if file.exists():
            mask = np.load(file, allow_pickle=False)
            audit = json.loads((out / f'mask_{graph}_audit.json').read_text())
            if repeated.digest(file) != audit['file_sha256']:
                raise ValueError('Saved mask checksum mismatch')
        else:
            mask, audit = shuffle_membership(projection, graph)
            mask = mask.astype(np.uint8)
            np.save(file, mask)
            audit['file_sha256'] = repeated.digest(file)
            repeated.save_json(out / f'mask_{graph}_audit.json', audit)
        if mask.shape != projection.shape or not np.isin(mask, [0, 1]).all():
            raise ValueError('Invalid mask')
        np.testing.assert_array_equal(mask.sum(0), (projection > 0).sum(0))
        np.testing.assert_array_equal(mask.sum(1), (projection > 0).sum(1))
        masks[graph] = mask
    records, baselines = [], []
    torch.set_num_threads(a.threads)
    with threadpool_limits(limits=a.threads):
        for split_seed in a.split_seeds:
            folder = out / f'split_{split_seed}'
            folder.mkdir(exist_ok=True)
            split = repeated.compound_split(data.names, split_seed)
            split.to_csv(folder / 'compound_split.csv', index=False)
            def rows(stage):
                return np.flatnonzero(data.pairs.compound.isin(split.loc[split.stage.eq(stage), 'compound']))
            itr, validation, testing = rows('inner_train'), rows('inner_validation'), rows('outer_test')
            training = np.sort(np.r_[itr, validation])
            _, baseline, detail = repeated.training_mean_response(data, training, testing)
            baseline['split_seed'] = split_seed
            baselines.append(baseline)
            detail.to_csv(folder / 'baseline_response_conditions.csv', index=False)
            pd.DataFrame(baselines).to_csv(out / 'response_baselines.csv', index=False)
            inner = prepare_stage(data, itr, validation, projection, pathways, provenance)
            outer = prepare_stage(data, training, testing, projection, pathways, provenance)
            joblib.dump(inner['preprocess'], folder / 'inner_preprocess.joblib')
            joblib.dump(outer['preprocess'], folder / 'refit_preprocess.joblib')
            for job in [j for j in planned if j['split_seed'] == split_seed]:
                name = f"{job['family']}_{job['pathway_mode']}_model{job['model_seed']}_graph{job['graph_seed']}"
                task = folder / name
                task.mkdir(exist_ok=True)
                marker = task / 'complete.json'
                if marker.exists():
                    record = json.loads(marker.read_text())
                    if any(record.get(k) != v for k, v in job.items()):
                        raise ValueError('Completed job identity differs')
                    for name in ['predictions.npz', 'model.pt', 'metrics.json', 'selection.json']:
                        if not (Path(record['attempt']) / name).is_file():
                            raise ValueError('Completed artifact missing')
                else:
                    attempt = task / ('attempt_' + uuid.uuid4().hex)
                    attempt.mkdir()
                    kind = {v: k for k, v in LABELS.items()}[job['family']]
                    mode = 'none' if job['pathway_mode'] == 'none' else 'sparse'
                    membership = masks[job['graph_seed']] if job['pathway_mode'] == 'random_sparse' else projection
                    print('Selecting on INNER validation:', job, flush=True)
                    choice = select_candidate(kind, data, inner, itr, validation, membership,
                                              mode, job['model_seed'], attempt)
                    repeated.save_json(attempt / 'selection.json', choice)
                    metrics = final_refit(kind, data, outer, training, testing, membership,
                                          mode, job['model_seed'], choice, attempt)
                    record = dict(**metrics, **job, structure=True, attempt=str(attempt))
                    repeated.save_json(attempt / 'metrics.json', record)
                    repeated.save_json(marker, record)
                records.append(record)
                pd.DataFrame(records).to_csv(out / 'results.csv', index=False)
    table = pd.DataFrame(records)
    if len(table) != len(planned):
        raise ValueError('Incomplete execution')
    per_unit, paired, split_paired, summary = repeated.summarize(table)
    for name, df in [('per_unit', per_unit), ('paired_differences', paired),
                     ('split_differences', split_paired), ('descriptive_summary', summary)]:
        df.to_csv(out / f'{name}.csv', index=False)
    against = per_unit.merge(pd.DataFrame(baselines), on='split_seed', validate='many_to_one')
    against['response_delta_vs_training_mean'] = against.response_rmse - against.training_mean_response_rmse
    against['response_delta_vs_zero'] = against.response_rmse - against.zero_response_rmse
    against.to_csv(out / 'models_vs_response_baselines.csv', index=False)
    repeated.save_json(out / 'completed.json', dict(completed_jobs=len(records)))


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['prepared', 'structure', 'annotation', 'out']:
        p.add_argument('--' + name, type=Path, required=True)
    p.add_argument('--families', nargs='+', choices=list(LABELS), default=['mlp'])
    p.add_argument('--split-seeds', nargs='+', type=int, default=[42, 43, 44, 45, 46])
    p.add_argument('--model-seeds', nargs='+', type=int, default=[42])
    p.add_argument('--graph-seeds', nargs='+', type=int, default=[20261102])
    p.add_argument('--threads', type=int, default=4)
    p.add_argument('--execute', action='store_true')
    return p


def main(argv=None):
    import repeated_pathway_validation as repeated
    p = parser()
    a = p.parse_args(argv)
    for field in ['families', 'split_seeds', 'model_seeds', 'graph_seeds']:
        values = getattr(a, field)
        if len(set(values)) != len(values):
            p.error('Duplicate ' + field)
    if a.threads < 1 or any(x < 0 or x >= 2**32 for field in
                           ['split_seeds', 'model_seeds', 'graph_seeds'] for x in getattr(a, field)):
        p.error('Invalid threads or seeds')
    planned = repeated.jobs(a.split_seeds, a.model_seeds, a.graph_seeds, [LABELS[f] for f in a.families])
    print(json.dumps(dict(mode='execute' if a.execute else 'plan_only', jobs=len(planned),
                         inner_candidate_fits=len(planned) * len(DEFAULT_PROTOCOL['learning_rates']),
                         final_refits=len(planned), fits_up_to=len(planned) * (len(DEFAULT_PROTOCOL['learning_rates']) + 1), protocol=DEFAULT_PROTOCOL), indent=2))
    if a.execute:
        execute(a, planned)
    return planned


if __name__ == '__main__':
    main()
