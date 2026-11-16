#!/usr/bin/env python3
"""Delta z profile (main) and per-structure netTE z profiles (supplementary).

    python make_profile_figures.py
"""
import os
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'data')
OUT  = os.path.join(HERE, '..', 'figures')

FS = dict(axis_label=17, tick=15, annot=14, legend=15, panel_letter=22, struct=16)
mpl.rcParams.update({
    'font.size': FS['tick'], 'axes.labelsize': FS['axis_label'],
    'xtick.labelsize': FS['tick'], 'ytick.labelsize': FS['tick'],
    'legend.fontsize': FS['legend'], 'axes.spines.top': False,
    'axes.spines.right': False, 'figure.facecolor': 'white',
    'savefig.facecolor': 'white', 'pdf.fonttype': 42, 'ps.fonttype': 42,
})
BLUE, RED, GREY, ORANGE = '#1f6fb4', '#b2182b', '#8a8a8a', '#d9731a'
DOMAINS = [(1, 390, 'TMD1'), (390, 650, 'NBD1'), (650, 860, 'R'),
           (860, 1170, 'TMD2'), (1170, 1480, 'NBD2')]

def domains(ax, ytop):
    for x0, x1, nm in DOMAINS:
        if nm in ('NBD1', 'NBD2'):
            ax.axvspan(x0, x1, color='#f4f4f4', zorder=0)
        ax.axvline(x0, color='0.86', lw=0.9, zorder=0)
        ax.text((x0 + x1) / 2, ytop, nm, ha='center', va='top',
                fontsize=FS['annot'], color='0.35')

# ======================= MAIN: Delta z profile ==============================
DZ  = pd.read_csv(f'{DATA}/deltaz_consensus.csv', index_col=0)['delta_z']
DZ.index = DZ.index.astype(int)
LAB = pd.read_csv(f'{DATA}/key_variant_labels.csv')

fig, ax = plt.subplots(figsize=(15.5, 6.4))
ytop = DZ.max() + 2.8
domains(ax, ytop)
ax.axhline(0, color='black', lw=1.0, zorder=1)
ax.fill_between(DZ.index, 0, DZ.values, where=DZ.values > 0, color=BLUE, alpha=0.30,
                interpolate=True, zorder=2)
ax.fill_between(DZ.index, 0, DZ.values, where=DZ.values < 0, color=ORANGE, alpha=0.30,
                interpolate=True, zorder=2)
ax.plot(DZ.index, DZ.values, lw=0.9, color='0.25', zorder=3)
for k, (_, r) in enumerate(LAB.iterrows()):
    res = int(r['residue'])
    if res not in DZ.index:
        continue
    y = DZ.loc[res]
    ax.scatter([res], [y], s=70, color=RED, zorder=5, edgecolors='white', linewidths=1.1)
    step = k % 3
    dy = ([34, 58, 82][step]) if y > 0 else (-[38, 62, 86][step])
    ax.annotate(r['variant'], xy=(res, y), xytext=(0, dy), textcoords='offset points',
                ha='center', fontsize=FS['annot'], color=RED,
                arrowprops=dict(arrowstyle='-', color=RED, lw=0.9, alpha=0.75))
ax.set_xlim(1, 1480); ax.set_ylim(DZ.min() - 1.4, ytop + 0.6)
ax.set_xlabel('CFTR residue number')
ax.set_ylabel('$\\Delta$z  (NBD-dimerized $-$ NBD-separated)')
ax.text(755, ytop - 1.3, 'positive: more sink-like when dimerized', ha='center',
        fontsize=FS['annot'], color=BLUE, va='top')
ax.text(1010, DZ.min() - 0.9, 'negative: more source-like when dimerized', ha='center',
        fontsize=FS['annot'], color=ORANGE, va='bottom')
fig.tight_layout()
fig.savefig(f'{OUT}/figure_deltaz_profile.png', dpi=300, bbox_inches='tight')
fig.savefig(f'{OUT}/figure_deltaz_profile.pdf', bbox_inches='tight')
plt.close(fig)

# ============ SUPPLEMENTARY: netTE z profile, every structure ===============
Z = pd.read_csv(f'{DATA}/netTE_z_per_structure.csv', index_col=0)
Z.index = Z.index.astype(int)
CF = pd.read_csv(f'{DATA}/cf_causing_residues.csv')
cfres = set(CF['res'])
ORDER = ['9mxl', '9dw9', '8eio', '6msm', '6o2p', '9dw8', '8v81', '5uak', '9dw5', '9dw7']
META = {'9mxl': '2.1 Å, NBD-dimerized', '9dw9': '2.8 Å, NBD-dimerized',
        '8eio': '2.8 Å, NBD-dimerized', '6msm': '3.2 Å, NBD-dimerized',
        '6o2p': '3.3 Å, NBD-dimerized', '9dw8': '3.5 Å, NBD-dimerized',
        '8v81': '3.6 Å, NBD-dimerized', '5uak': '3.87 Å, NBD-separated',
        '9dw5': '3.8 Å, NBD-separated', '9dw7': '6.0 Å, NBD-separated'}

fig, axes = plt.subplots(len(ORDER), 1, figsize=(15.5, 2.9 * len(ORDER)), sharex=True)
for ax, col in zip(axes, ORDER):
    z = Z[col].dropna()
    ytop = z.max() + 1.6
    domains(ax, ytop)
    ax.axhline(0, color=GREY, lw=0.9, ls=':', zorder=1)
    ax.plot(z.index, z.values, lw=1.0, color=BLUE, zorder=2)
    hit = [r for r in z.index if r in cfres]
    ax.scatter(hit, z.loc[hit], s=30, color=RED, zorder=4,
               edgecolors='white', linewidths=0.6, label='CF-causing residue')
    ax.set_xlim(1, 1480); ax.set_ylim(z.min() - 0.5, ytop + 0.5)
    ax.set_ylabel('netTE (z)')
    # structure name on its own line above the axes - never overlaps the trace
    ax.set_title('%s — %s' % (col.upper(), META[col]), loc='left',
                 fontsize=FS['struct'], fontweight='bold', pad=8)
axes[0].legend(frameon=False, loc='upper center', bbox_to_anchor=(0.52, 1.02))
axes[-1].set_xlabel('CFTR residue number')
fig.tight_layout(h_pad=2.2)
fig.savefig(f'{OUT}/figureS_profiles_all_structures.png', dpi=300, bbox_inches='tight')
fig.savefig(f'{OUT}/figureS_profiles_all_structures.pdf', bbox_inches='tight')
plt.close(fig)
print('wrote figures to', os.path.normpath(OUT))
