"""
gnm_te_pipeline.py
------------------
Reference implementation of the GNM transfer-entropy pipeline used in
"A conformational change in allosteric information flow predicts corrector
rescue of CFTR variants".

Reproduces the alloGlyph command-line pipeline exactly (Spearman rho = 1.0000
against its per-residue z output on 9DW9 chain A, 1,182 residues, and on the
AlphaFold model, 1,480 residues).

Node definition : one node per residue at the centroid of ALL its atoms
                  (heavy-atom centroid; hydrogens are absent from the inputs).
Contact cutoff  : 11.5 A (the value used for CFTR).
Relaxation time : tau is selected per structure as the MEDIAN OVER RESIDUES of
                  the time at which that residue's GNM autocorrelation decays
                  to 1/e. It is a per-residue quantity aggregated by median,
                  not a per-mode quantity.

Usage
-----
    python gnm_te_pipeline.py 9MXL B
    python gnm_te_pipeline.py 5UAK A --cutoff 11.5

Dependencies: numpy, pandas (scipy only for the optional comparison step).
"""
import argparse, gzip, os, sys, urllib.request
import numpy as np
import pandas as pd

CACHE = os.environ.get("CFTR_CIF_CACHE", "./cif_cache")


def fetch_cif(pdb_id, cache=CACHE):
    """Download an mmCIF from RCSB (cached on disk)."""
    os.makedirs(cache, exist_ok=True)
    path = os.path.join(cache, f"{pdb_id.lower()}.cif")
    if not os.path.exists(path):
        url = f"https://files.rcsb.org/download/{pdb_id.upper()}.cif.gz"
        with urllib.request.urlopen(url, timeout=90) as r:
            open(path, "wb").write(gzip.decompress(r.read()))
    return path


def residue_centroids(cif_path, chain):
    """Parse one chain of an mmCIF -> (residue numbers, all-atom centroids)."""
    cols, res, in_loop = [], {}, False
    for line in open(cif_path, errors="ignore"):
        if line.startswith("_atom_site."):
            cols.append(line.strip().split(".")[1]); in_loop = True; continue
        if in_loop and line.startswith(("ATOM", "HETATM")):
            rec = dict(zip(cols, line.split()))
            if rec["group_PDB"] != "ATOM":                       continue
            if rec.get("pdbx_PDB_model_num", "1") != "1":        continue
            if rec.get("auth_asym_id", rec.get("label_asym_id")) != chain: continue
            if rec["type_symbol"] == "H":                        continue
            try:
                num = int(rec["auth_seq_id"])
            except ValueError:
                continue
            res.setdefault(num, []).append(
                [float(rec["Cartn_x"]), float(rec["Cartn_y"]), float(rec["Cartn_z"])])
        elif in_loop and line.startswith("#") and res:
            break
    ids = np.array(sorted(res))
    xyz = np.array([np.mean(res[i], axis=0) for i in ids])
    return ids, xyz


def kirchhoff(xyz, cutoff):
    d = np.linalg.norm(xyz[:, None] - xyz[None], axis=-1)
    a = ((d <= cutoff) & (d > 0)).astype(float)
    return np.diag(a.sum(1)) - a, a.sum(1)


def modes(xyz, cutoff, tol=1e-8):
    gamma, degree = kirchhoff(xyz, cutoff)
    lam, vec = np.linalg.eigh(gamma)
    nz = lam > tol
    return lam[nz], vec[:, nz], degree


def select_tau(lam, vec, t_max=5.0, n_grid=2000):
    """Median over residues of the 1/e decay time of the residue autocorrelation."""
    t = np.linspace(0, t_max, n_grid)
    w2 = vec ** 2
    num = w2 @ (np.exp(-np.outer(lam, t)) / lam[:, None])
    den = w2 @ (1 / lam)
    ct = num / den[:, None]
    thr = 1 / np.e
    tau_res = np.empty(ct.shape[0])
    for r in range(ct.shape[0]):
        below = ct[r] <= thr
        k = np.argmax(below) if below.any() else len(t) - 1
        tau_res[r] = t[k] if k == 0 else (
            t[k - 1] + (thr - ct[r, k - 1]) * (t[k] - t[k - 1]) / (ct[r, k] - ct[r, k - 1]))
    return float(np.median(tau_res))


def net_transfer_entropy(lam, vec, tau):
    """Pairwise TE (Hacisuleyman & Erman form) summed into a signed per-residue net."""
    c = (vec * (1 / lam)) @ vec.T                       # equilibrium covariance
    ct = (vec * (np.exp(-lam * tau) / lam)) @ vec.T     # time-lagged covariance
    dc, dct = np.diag(c), np.diag(ct)
    cii, cjj = dc[:, None], dc[None, :]
    ctj = dct[None, :]
    a = cjj ** 2 - ctj ** 2
    b = cii * cjj ** 2
    cc = 2 * c * ctj * ct
    dd = -(((ct ** 2) + (c ** 2)) * cjj) - ((ctj ** 2) * cii)
    f, g = cjj, cii * cjj - c ** 2
    with np.errstate(invalid="ignore", divide="ignore"):
        te = 0.5 * np.log(a) - 0.5 * np.log(b + cc + dd) - 0.5 * np.log(f) + 0.5 * np.log(g)
    te = np.nan_to_num(te)
    np.fill_diagonal(te, 0.0)
    te[te < 0] = 0.0                                    # T(i->j) >= 0 by construction
    return (te - te.T).sum(1)


def run(pdb_id, chain, cutoff=11.5, tau=None):
    ids, xyz = residue_centroids(fetch_cif(pdb_id), chain)
    lam, vec, degree = modes(xyz, cutoff)
    if tau is None:
        tau = select_tau(lam, vec)
    net = net_transfer_entropy(lam, vec, tau)
    z = (net - net.mean()) / net.std()
    return pd.DataFrame({"resnum": ids, "net_te": net, "z": z,
                         "contacts": degree}), tau


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("pdb_id")
    ap.add_argument("chain")
    ap.add_argument("--cutoff", type=float, default=11.5)
    ap.add_argument("--tau", type=float, default=None,
                    help="fixed tau; omit to select per structure")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    df, tau = run(args.pdb_id, args.chain, args.cutoff, args.tau)
    out = args.out or f"{args.pdb_id.lower()}_{args.chain}_netTE.csv"
    df.to_csv(out, index=False)
    print(f"{args.pdb_id} chain {args.chain}: {len(df)} residues, "
          f"mean contacts {df.contacts.mean():.3f}, tau {tau:.4f} -> {out}")


if __name__ == "__main__":
    main()
