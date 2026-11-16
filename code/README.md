# Code and data — CFTR Δz analysis

Code and data for *"A conformational change in allosteric information
flow predicts corrector rescue of CFTR variants"* (Aysima Hacisuleyman, Department of Computational Biology,
University of Lausanne).

## Layout

    code/      analysis and figure scripts (this folder)
    ../data/   processed tables; every figure panel reads from here
    ../figures/ rendered figures (PNG + PDF, 300 dpi, Type-42 fonts)

## Scripts

| file | produces |
|---|---|
| `gnm_te_pipeline.py` | reference implementation: fetches a structure, builds the Gaussian network, selects τ, returns per-residue net transfer entropy and its z-score |
| `make_figures.py` | Figures 1, 3, 4 and 5 |
| `make_profile_figures.py` | Figure S1 (per-structure netTE profiles) and the standalone Δz profile |
| `render_deltaz_structure.pml` | Figure 2 (PyMOL; reads `../data/cftr_9mxl_deltaz.pdb`) |
| `netTE_substitution_score.py` | substitution-level scoring used in the benchmark (Figure 3d) |

## Reproducing the core quantity

    python gnm_te_pipeline.py 9MXL B

prints `1173 residues, mean contacts 26.121, tau 0.0581`, the values reported in
Methods, and writes the per-residue table.

The implementation reproduces the alloGlyph command-line pipeline exactly:
Spearman ρ = 1.0000 against its per-residue z output on 9DW9 chain A (1,182
residues) and on the AlphaFold model AF-P13569-F1-model_v6 (1,480 residues).

### Parameters

* **Nodes** — one per residue at the centroid of all its atoms (heavy-atom
  centroid; inputs carry no hydrogens).
* **Contact cutoff** — 11.5 Å.
* **τ** — selected per structure as the *median over residues* of the time at
  which that residue's GNM autocorrelation decays to 1/e. A per-residue quantity
  aggregated by median, not a per-mode quantity. Values: 0.0554–0.0585 across the
  seven NBD-dimerized structures, 0.0971–0.1016 across the three NBD-separated
  structures, 0.0697 for the AlphaFold model.
* **Sign** — netTE_i = Σ_j (T_i→j − T_j→i); positive for sources, negative for
  sinks. Pairwise T(i→j) ≥ 0 by construction.

## Figure 2 (structure)

`../data/cftr_9mxl_deltaz.pdb` is 9MXL chain B with Δz written into the B-factor
column (sentinel 99.00 for the 38 residues that have no Δz because they are
unresolved in the NBD-separated structures). It opens in any viewer and colours
by Δz directly. Render with:

    pymol -cq render_deltaz_structure.pml

## Data files

`../data/figNx_*.csv` are named for the figure panel they feed, e.g.
`fig4d_forest.csv` → Figure 4d. Other tables:

| file | contents |
|---|---|
| `netTE_z_per_structure.csv` | per-residue netTE z for all ten experimental structures (1,103 residues resolved in all) |
| `deltaz_consensus.csv` | Δz for the 1,135 residues resolved in both states |
| `key_variants_deltaz.csv` | clinically important variants with z, Δz, percentile, CFTR2 class, corrector gain |
| `deltaz_matched_tau.csv` | Δz–rescue correlations recomputed at fixed τ |
| `tau_sensitivity.csv`, `source_sink_inversion.csv` | τ dependence of the sink/source assignment |
| `enrichment_dimerised_consensus.csv` | tier enrichment by variant class |
| `clinvar_mapping_validation.csv` | genomic→residue mapping checks |

## Dependencies

Python 3.11 with numpy, scipy, pandas, matplotlib, statsmodels.
Figure 2 additionally requires PyMOL (tested with pymol-open-source 3.x).

## External data

Structures from the RCSB PDB (9MXL, 9DW9, 9DW8, 9DW7, 9DW5, 8EIO, 8V81, 6MSM,
6O2P, 5UAK) and the AlphaFold Protein Structure Database (AF-P13569-F1-model_v6);
fetched at run time. Variant data from CFTR2 (release 30 January 2026), ClinVar
via the gnomAD r4 mirror (snapshot 11 October 2026) and AlphaMissense. Functional
measurements are from the four published datasets cited in the Methods and are
not redistributed here.
