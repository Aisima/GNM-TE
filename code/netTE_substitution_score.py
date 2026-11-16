#!/usr/bin/env python3
"""
Per-variant score = per-residue netTE + burial + substitution chemistry.

Why: netTE gives ONE value per residue, but pathogenicity differs between
substitutions at the same position (R117H vs R117P). A per-residue score is
therefore capped -- on CFTR a perfect per-residue oracle reaches AUC 0.910.
Adding four amino-acid features lifts cross-validated AUC from 0.747 to 0.855,
and netTE still contributes +0.016 AUC over burial+chemistry alone
(95% CI 0.010-0.022, residue-level bootstrap).

Usage
-----
    python netTE_substitution_score.py track.csv variants.csv out.csv
       [--pdb structure.cif]      # for the burial term; else burial is dropped

track.csv     alloGlyph track (needs: chain, resnum, z)
variants.csv  one row per variant: resnum, wt, mut   (1-letter amino acids)

Refit on your own labels (recommended -- the shipped coefficients were fitted
on CFTR AlphaMissense calls and are a starting point, not a trained predictor):
    from netTE_substitution_score import featurize, fit
    coef = fit(featurized_df, labels)
"""
import argparse, sys
import numpy as np
import pandas as pd

KD = dict(A=1.8, R=-4.5, N=-3.5, D=-3.5, C=2.5, Q=-3.5, E=-3.5, G=-0.4, H=-3.2,
          I=4.5, L=3.8, K=-3.9, M=1.9, F=2.8, P=-1.6, S=-0.8, T=-0.7, W=-0.9,
          Y=-1.3, V=4.2)                                   # Kyte-Doolittle
VOL = dict(A=88.6, R=173.4, N=114.1, D=111.1, C=108.5, Q=143.8, E=138.4, G=60.1,
           H=153.2, I=166.7, L=166.7, K=168.6, M=162.9, F=189.9, P=112.7,
           S=89.0, T=116.1, W=227.8, Y=193.6, V=140.0)     # A^3
CHG = dict(D=-1, E=-1, K=1, R=1, H=0.5)

# Fitted on 18,914 CFTR substitutions (AlphaMissense labels, 11.5 A consensus
# netTE over 7 NBD-dimerised structures). Refit for other systems.
COEF = dict(const=-4.7746, z=0.6902, cn=0.1835, dhyd=0.3191,
            dvol=0.0233, dchg=0.8990, pg=0.7406)


def burial_from_structure(path, chain, cutoff=10.0):
    """CA neighbour count within `cutoff` A. Accepts .cif/.pdb (plain or .gz)."""
    import gzip
    op = gzip.open if str(path).endswith('.gz') else open
    xyz, res = [], []
    with op(path, 'rt') as fh:
        lines = fh.read().splitlines()
    if any(l.startswith('_atom_site.') for l in lines):           # mmCIF
        hdr = [l.split('.')[1].strip() for l in lines if l.startswith('_atom_site.')]
        ix = {k: hdr.index(k) for k in ('label_atom_id', 'auth_asym_id',
                                        'auth_seq_id', 'Cartn_x', 'Cartn_y', 'Cartn_z')}
        for l in lines:
            if not l.startswith('ATOM'):
                continue
            f = l.split()
            if f[ix['label_atom_id']].strip('"') != 'CA' or f[ix['auth_asym_id']] != chain:
                continue
            res.append(int(f[ix['auth_seq_id']]))
            xyz.append([float(f[ix[c]]) for c in ('Cartn_x', 'Cartn_y', 'Cartn_z')])
    else:                                                          # PDB
        for l in lines:
            if l.startswith('ATOM') and l[12:16].strip() == 'CA' and l[21] == chain:
                res.append(int(l[22:26]))
                xyz.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
    if not res:
        raise SystemExit(f'no CA atoms for chain {chain} in {path}')
    P = np.asarray(xyz)
    D = np.linalg.norm(P[:, None] - P[None], axis=-1)
    return pd.Series(((D < cutoff) & (D > 0)).sum(1).astype(float), index=res, name='cn')


def featurize(variants, track, burial=None, chain=None):
    """variants: resnum/wt/mut. track: alloGlyph track. -> feature frame."""
    t = track.copy()
    if chain is not None and 'chain' in t.columns:
        t = t[t.chain == chain]
    zmap = t.set_index('resnum')['z']
    v = variants.copy()
    bad = ~v.wt.isin(KD) | ~v.mut.isin(KD)
    if bad.any():
        print(f'[warn] dropping {int(bad.sum())} rows with non-standard residues', file=sys.stderr)
        v = v[~bad]
    v['z'] = -v.resnum.map(zmap)                 # sign: sink-like = higher risk
    v['cn'] = v.resnum.map(burial) if burial is not None else 0.0
    v['dhyd'] = (v.wt.map(KD) - v.mut.map(KD)).abs()
    v['dvol'] = (v.wt.map(VOL) - v.mut.map(VOL)).abs()
    v['dchg'] = (v.wt.map(CHG).fillna(0) - v.mut.map(CHG).fillna(0)).abs()
    v['pg'] = (v.wt.isin(['P', 'G']) | v.mut.isin(['P', 'G'])).astype(int)
    return v


def score(feat, coef=COEF, use_burial=True):
    terms = ['z', 'cn', 'dhyd', 'dvol', 'dchg', 'pg'] if use_burial else \
            ['z', 'dhyd', 'dvol', 'dchg', 'pg']
    lin = coef['const'] + sum(coef[k] * feat[k].astype(float) for k in terms)
    return 1.0 / (1.0 + np.exp(-lin))


def fit(feat, labels, use_burial=True):
    """Refit the logistic model on your own 0/1 labels. Returns a coef dict."""
    import statsmodels.api as sm
    terms = ['z', 'cn', 'dhyd', 'dvol', 'dchg', 'pg'] if use_burial else \
            ['z', 'dhyd', 'dvol', 'dchg', 'pg']
    X = sm.add_constant(feat[terms].astype(float))
    m = sm.Logit(labels, X).fit(disp=0)
    return {k: float(val) for k, val in m.params.items()}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('track'); ap.add_argument('variants'); ap.add_argument('out')
    ap.add_argument('--pdb', help='structure for the burial term (.cif/.pdb, .gz ok)')
    ap.add_argument('--chain', default=None, help='chain id; default = first in track')
    a = ap.parse_args()

    track = pd.read_csv(a.track)
    variants = pd.read_csv(a.variants)
    chain = a.chain or (track.chain.iloc[0] if 'chain' in track.columns else None)

    burial = burial_from_structure(a.pdb, chain) if a.pdb else None
    if burial is None:
        print('[warn] no --pdb given: scoring without the burial term', file=sys.stderr)

    feat = featurize(variants, track, burial, chain)
    feat['score'] = score(feat, use_burial=burial is not None)
    feat.sort_values('score', ascending=False).to_csv(a.out, index=False)
    print(f'wrote {a.out}: {len(feat)} variants, score range '
          f'{feat.score.min():.3f}-{feat.score.max():.3f}')


if __name__ == '__main__':
    main()
