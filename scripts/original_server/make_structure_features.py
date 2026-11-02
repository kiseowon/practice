"""Build 2D / 3D structure feature tables from a SMILES CSV (needs rdkit: pip install rdkit).
Usage: python make_structure_features.py --smiles compound_pubchem_smiles.csv --outdir struct_feats [--name_col compound --smiles_col smiles]
Outputs (one row per compound, first column = compound):
  struct_morgan2048.csv  (2D, ECFP4 bits)   struct_maccs.csv (2D, 167 keys)   struct_rdkit2d.csv (2D, ~210 RDKit descriptors)
  struct_3d.csv (3D: ETKDG+MMFF conformer -> PMI/NPR/Rg/asphericity/... + RDF + MORSE + WHIM + AUTOCORR3D)
"""
import argparse, os, re, numpy as np, pandas as pd, warnings
warnings.filterwarnings("ignore")
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, MACCSkeys, Descriptors, rdMolDescriptors as rd
from rdkit.Chem.MolStandardize import rdMolStandardize
RDLogger.DisableLog("rdApp.*")
ap = argparse.ArgumentParser(); ap.add_argument("--smiles", required=True); ap.add_argument("--outdir", default="struct_feats")
ap.add_argument("--name_col"); ap.add_argument("--smiles_col"); ap.add_argument("--seeds", type=int, default=5); a = ap.parse_args()
os.makedirs(a.outdir, exist_ok=True)
d = pd.read_csv(a.smiles, low_memory=False); print("columns:", list(d.columns))
sc = a.smiles_col or next((c for c in d.columns if re.search("smiles", c, re.I)), None)
nc = a.name_col or next((c for c in d.columns if re.search("compound|name", c, re.I) and c != sc), d.columns[0])
print("name col:", nc, "| smiles col:", sc); assert sc, "no SMILES column found; pass --smiles_col"
chooser = rdMolStandardize.LargestFragmentChooser()
def prep(smi):
    m = Chem.MolFromSmiles(str(smi))
    return chooser.choose(m) if m is not None else None
mols = [prep(s) for s in d[sc]]; ok = [m is not None for m in mols]; print("parsed:", sum(ok), "of", len(ok))

def table(feats, name):
    df = pd.DataFrame(feats); df.insert(0, "compound", d[nc].values)
    X = df.drop(columns="compound").replace([np.inf, -np.inf], np.nan)
    X = X.fillna(X.median()); X = X.loc[:, X.std() > 0]; df = pd.concat([df[["compound"]], X], axis=1)
    df.to_csv(os.path.join(a.outdir, name), index=False); print(name, "->", df.shape)

n = len(mols); z = lambda k: [np.nan] * k
table([list(AllChem.GetMorganFingerprintAsBitVect(m, 2, 2048)) if m else z(2048) for m in mols], "struct_morgan2048.csv")
table([list(MACCSkeys.GenMACCSKeys(m)) if m else z(167) for m in mols], "struct_maccs.csv")
names = [x[0] for x in Descriptors._descList]
table([[f(m) for _, f in Descriptors._descList] if m else z(len(names)) for m in mols], "struct_rdkit2d.csv")

def d3(m):
    """3D descriptors averaged over several ETKDG conformers (lowest-energy MMFF conformer per seed)"""
    if m is None: return None
    mh = Chem.AddHs(m); vs = []
    for seed in range(a.seeds):
        c = Chem.Mol(mh)
        ps = AllChem.ETKDGv3(); ps.randomSeed = seed + 1
        if AllChem.EmbedMolecule(c, ps) != 0: continue
        try: AllChem.MMFFOptimizeMolecule(c, maxIters=500)
        except Exception: pass
        f = [rd.CalcPMI1(c), rd.CalcPMI2(c), rd.CalcPMI3(c), rd.CalcNPR1(c), rd.CalcNPR2(c), rd.CalcRadiusOfGyration(c),
             rd.CalcInertialShapeFactor(c), rd.CalcEccentricity(c), rd.CalcAsphericity(c), rd.CalcSpherocityIndex(c), rd.CalcPBF(c)]
        f += list(rd.CalcRDF(c)) + list(rd.CalcMORSE(c)) + list(rd.CalcWHIM(c)) + list(rd.CalcAUTOCORR3D(c))
        vs.append(f)
    return np.mean(vs, 0) if vs else None
rows = [d3(m) for m in mols]; k = next(len(r) for r in rows if r is not None)
table([r if r is not None else z(k) for r in rows], "struct_3d.csv")
print("3D failed for", sum(r is None for r in rows), "compounds")
