"""Reconstructed preparation utilities, NOT the original server mapping pipeline.

Never establishes provenance of the existing final inputs. Uses fresh output dirs.
No training. External retrieval occurs only when its subcommand is invoked.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import time
from urllib.request import urlopen


def write_csv(path, fields, rows):
    with path.open('w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def read_csv(path):
    with path.open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))


def smiles(a):
    import pubchempy as pcp
    names = [r['compound'].strip() for r in read_csv(a.compounds)]
    if len(names) != len(set(names)):
        raise ValueError('Duplicate compound names')
    rows = []
    for name in names:
        row = dict(compound=name, pubchem_cid='', smiles='', status='NOT_FOUND')
        if name in {'acetamidofluorene', 'naphthyl isothiocyanate'}:
            # Original manual CID/SMILES values were not supplied; never invent them.
            row['status'] = 'MANUAL_REQUIRED'
        else:
            for attempt in range(4):  # initial attempt + at most 3 retries
                try:
                    found = pcp.get_compounds(name, 'name')
                    if found:
                        c = found[0]
                        row.update(pubchem_cid=c.cid, smiles=c.connectivity_smiles,
                                   status='FIRST_HIT_UNVERIFIED')
                    break
                except Exception as exc:
                    row['status'] = 'ERROR_' + type(exc).__name__
                    if attempt < 3:
                        time.sleep(1)
        rows.append(row)
    write_csv(a.out / 'compound_pubchem_candidates.csv',
              ['compound', 'pubchem_cid', 'smiles', 'status'], rows)
    print('Candidates only. Review identity and add reviewed=true before Morgan generation.')


def morgan(a):
    from rdkit import Chem, rdBase
    from rdkit.Chem import AllChem
    from rdkit.Chem.MolStandardize import rdMolStandardize
    chooser = rdMolStandardize.LargestFragmentChooser()
    rows, audit, seen = [], [], set()
    for r in read_csv(a.reviewed_smiles):
        name = r['compound'].strip()
        if name in seen:
            raise ValueError('Duplicate compound: ' + name)
        seen.add(name)
        if r.get('reviewed', '').strip().lower() != 'true':
            audit.append(dict(compound=name, status='UNREVIEWED_SKIPPED'))
            continue
        mol = Chem.MolFromSmiles(r['smiles'])
        if mol is None:
            audit.append(dict(compound=name, status='PARSE_FAILED'))
            continue
        mol = chooser.choose(mol)
        bits = AllChem.GetMorganFingerprintAsBitVect(mol, a.radius, nBits=2048)
        rows.append(dict(compound=name, **{f'morgan_{i}': int(bits[i]) for i in range(2048)}))
        audit.append(dict(compound=name, status='OK_RECONSTRUCTED'))
    write_csv(a.out / 'morgan_binary.csv', ['compound'] + [f'morgan_{i}' for i in range(2048)], rows)
    write_csv(a.out / 'structure_audit.csv', ['compound', 'status'], audit)
    return {'radius': a.radius, 'n_bits': 2048, 'rdkit': rdBase.rdkitVersion,
            'policy': 'reviewed only; largest fragment; no median imputation or constant-column removal'}


def kegg(a):
    import numpy as np
    # Current KEGG responses are saved verbatim. This cannot recreate an old snapshot.
    texts = {}
    for key, route in [('genes', 'list/rno'), ('pathways', 'list/pathway/rno'),
                       ('membership', 'link/pathway/rno')]:
        with urlopen('https://rest.kegg.jp/' + route, timeout=60) as response:
            raw = response.read()
        (a.out / f'kegg_{key}.txt').write_bytes(raw)
        texts[key] = raw.decode('utf-8')
    symbol_genes = {}
    for line in texts['genes'].splitlines():
        parts = line.split('\t')
        if len(parts) < 2:
            continue
        # KEGG gene list description begins with comma-separated aliases before ';'.
        for symbol in parts[-1].split(';')[0].split(','):
            symbol_genes.setdefault(symbol.strip(), set()).add(parts[0])
    mapping = {}
    for r in read_csv(a.gene_map):
        mapping.setdefault(r['PROBEID'], set()).add(r['SYMBOL'].strip())
    probes = [r['probe'] for r in read_csv(a.probes)]
    if len(probes) != len(set(probes)):
        raise ValueError('Duplicate probes')
    probe_genes, audit = [], []
    for probe in probes:
        candidates = set()
        symbols = mapping.get(probe, set())
        for symbol in symbols:
            candidates.update(symbol_genes.get(symbol, set()))
        gene = next(iter(candidates)) if len(candidates) == 1 else None
        probe_genes.append(gene)
        audit.append(dict(probe=probe, symbols='|'.join(sorted(symbols)),
                          kegg_gene=gene or '', candidate_genes='|'.join(sorted(candidates)),
                          status='MATCHED' if gene else 'AMBIGUOUS' if candidates else 'UNMATCHED'))
    membership = {}
    for line in texts['membership'].splitlines():
        gene, path = line.split('\t')
        membership.setdefault(path.removeprefix('path:'), set()).add(gene)
    names = dict(line.split('\t', 1) for line in texts['pathways'].splitlines())
    columns, pathways = [], []
    for path, genes in sorted(membership.items()):
        matched = genes & set(g for g in probe_genes if g)
        if len(matched) < a.min_genes or len(genes) > a.max_genes:
            continue
        mask = np.array([g in genes if g else False for g in probe_genes], dtype=float)
        columns.append(mask / mask.sum())
        pathways.append(dict(pathway=path, name=names.get(path, names.get('path:' + path, '')),
                             total_kegg_genes=len(genes), matched_genes=len(matched),
                             matched_probes=int(mask.sum())))
    if not columns:
        raise ValueError('No retained pathways; inspect downloaded responses and mapping audit')
    np.save(a.out / 'projection.npy', np.column_stack(columns))
    write_csv(a.out / 'probe_mapping_audit.csv',
              ['probe', 'symbols', 'kegg_gene', 'candidate_genes', 'status'], audit)
    write_csv(a.out / 'pathway_audit.csv',
              ['pathway', 'name', 'total_kegg_genes', 'matched_genes', 'matched_probes'], pathways)
    return {'species': 'Rattus norvegicus', 'probe_order': probes,
            'mapping_policy': 'exact SYMBOL/KEGG alias; only unique KEGG candidate retained',
            'minimum_matched_genes': a.min_genes, 'maximum_total_pathway_genes': a.max_genes,
            'sources': ['https://rest.kegg.jp/list/rno', 'https://rest.kegg.jp/list/pathway/rno',
                        'https://rest.kegg.jp/link/pathway/rno']}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest='command', required=True)
    s = sub.add_parser('smiles'); s.add_argument('--compounds', type=Path, required=True)
    m = sub.add_parser('morgan'); m.add_argument('--reviewed-smiles', type=Path, required=True)
    m.add_argument('--radius', type=int, required=True, help='Explicit NEW setting, not verified historical radius')
    k = sub.add_parser('kegg'); k.add_argument('--gene-map', type=Path, required=True)
    k.add_argument('--probes', type=Path, required=True)
    k.add_argument('--min-genes', type=int, default=5); k.add_argument('--max-genes', type=int, default=500)
    for child in (s, m, k):
        child.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    if a.command == 'morgan' and a.radius < 0:
        p.error('radius must be nonnegative')
    if a.command == 'kegg' and (a.min_genes < 1 or a.max_genes < a.min_genes):
        p.error('invalid pathway size bounds')
    a.out.mkdir(parents=True, exist_ok=False)
    details = globals()[a.command](a) or {}
    inputs = {}
    for field in ('compounds', 'reviewed_smiles', 'gene_map', 'probes'):
        path = getattr(a, field, None)
        if path:
            inputs[field] = {'file': str(path), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}
    record = {'reconstructed': True, 'historical_final_input_equivalence': 'UNVERIFIED',
              'command': a.command, 'input_files': inputs, 'created_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()), **details}
    (a.out / 'provenance.json').write_text(json.dumps(record, ensure_ascii=False, indent=2) + '\n')


if __name__ == '__main__':
    main()
