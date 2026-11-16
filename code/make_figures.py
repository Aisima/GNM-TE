#!/usr/bin/env python3
"""Regenerate Figures 2 and 3 for the CFTR netTE / Delta z manuscript.

Reads the per-panel source tables in ../data/ and writes PNGs to ../figures/.
No descriptive titles on any panel; panel letters only. Font sizes are set by FS.

    python make_figures.py
"""
import os
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'data')
OUT  = os.path.join(HERE, '..', 'figures')
os.makedirs(OUT, exist_ok=True)

# ---- typography ------------------------------------------------------------
FS = dict(axis_label=15, tick=13, annot=13, legend=13, panel_letter=20,
          heat_cell=13, ytick_forest=13)
mpl.rcParams.update({
    'font.size': FS['tick'], 'axes.labelsize': FS['axis_label'],
    'xtick.labelsize': FS['tick'], 'ytick.labelsize': FS['tick'],
    'legend.fontsize': FS['legend'], 'axes.titlesize': FS['panel_letter'],
    'axes.spines.top': False, 'axes.spines.right': False,
    'figure.facecolor': 'white', 'savefig.facecolor': 'white',
    'pdf.fonttype': 42, 'ps.fonttype': 42,
})
BLUE, ORANGE, GREY = '#1f6fb4', '#d9731a', '#8a8a8a'
RED, GREEN = '#b2182b', '#1b7837'

def panel(ax, letter):
    ax.text(-0.02, 1.06, letter, transform=ax.transAxes,
            fontsize=FS['panel_letter'], fontweight='bold', va='bottom', ha='right')

# ============================ FIGURE 1 =======================================
Cm  = pd.read_csv(f'{DATA}/fig1a_correlation_matrix.csv', index_col=0)
PR  = pd.read_csv(f'{DATA}/fig1b_consensus_profiles.csv')
GEO = pd.read_csv(f'{DATA}/fig1c_geometry.csv')

fig = plt.figure(figsize=(16.0, 6.4))
gs = fig.add_gridspec(1, 3, width_ratios=[0.95, 2.15, 0.58], wspace=0.32,
                      left=0.05, right=0.975)

# (a) cross-structure correlation matrix
ax = fig.add_subplot(gs[0, 0])
im = ax.imshow(Cm.values, cmap='viridis', vmin=0, vmax=1)
names = [c.upper() if c != 'AF' else 'AF' for c in Cm.columns]
ax.set_xticks(range(len(names))); ax.set_xticklabels(names, rotation=90, fontsize=FS['tick']-2)
ax.set_yticks(range(len(names))); ax.set_yticklabels(names, fontsize=FS['tick']-2)
nd = int((GEO['state'] == 'dimerized').sum()); nsep = int((GEO['state'] == 'separated').sum())
for b in (nd-0.5, nd+nsep-0.5):
    ax.axhline(b, color='white', lw=2.5); ax.axvline(b, color='white', lw=2.5)
cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
cb.set_label('Spearman $\\rho$', fontsize=FS['axis_label']); cb.ax.tick_params(labelsize=FS['tick']-1)
panel(ax, 'a')

# (b) consensus netTE profiles along the sequence
ax = fig.add_subplot(gs[0, 1])
ax.axhline(0, color=GREY, lw=1.0, ls=':')
ax.plot(PR['residue'], PR['dimerized_consensus'], lw=1.2, color=BLUE, label='NBD-dimerized')
ax.plot(PR['residue'], PR['separated_consensus'], lw=1.2, color=ORANGE, alpha=0.85,
        label='NBD-separated')
for x0, x1, nm in [(1, 390, 'TMD1'), (390, 650, 'NBD1'), (650, 860, 'R'),
                   (860, 1170, 'TMD2'), (1170, 1480, 'NBD2')]:
    ax.axvline(x0, color='0.85', lw=0.9, zorder=0)
    ax.text((x0+x1)/2, ax.get_ylim()[1]*0.97, nm, ha='center', va='top',
            fontsize=FS['annot']-1, color='0.35')
ax.set_xlabel('CFTR residue number')
ax.set_ylabel('netTE (z, consensus)')
ax.set_xlim(1, 1480)
ax.legend(frameon=False, loc='lower right', ncol=2)
panel(ax, 'b')

# (c) ATP-site separation by state
ax = fig.add_subplot(gs[0, 2])
cmap_state = {'dimerized': BLUE, 'separated': ORANGE, 'AlphaFold': GREEN}
xj = {'dimerized': 0, 'separated': 1, 'AlphaFold': 2}
rng = np.random.default_rng(7)
for _, r in GEO.iterrows():
    x = xj[r['state']] + rng.uniform(-0.13, 0.13)
    ax.scatter(x, r['atp_site_separation_A'], s=130, color=cmap_state[r['state']],
               edgecolors='white', linewidths=1.3, zorder=3)
ax.set_xticks([0, 1, 2])
ax.set_xticklabels(['Dimerized', 'Separated', 'AlphaFold'], fontsize=FS['tick']-1,
                   rotation=90, va='top')
ax.set_ylabel('ATP-site C$\\alpha$ separation (\u00c5)')
ax.set_xlim(-0.5, 2.5); ax.set_ylim(0, 33)
ax.axhspan(13.2, 27.5, color='#f2f2f2', zorder=0)
panel(ax, 'c')

fig.savefig(f'{OUT}/figure1_two_architectures.png', dpi=300, bbox_inches='tight')
fig.savefig(f'{OUT}/figure1_two_architectures.pdf', bbox_inches='tight')
plt.close(fig)

# ============================ FIGURE 2 (scope) ===============================
EN  = pd.read_csv(f'{DATA}/fig3a_enrichment.csv')
SC  = pd.read_csv(f'{DATA}/fig3b_scatter.csv')
TAU = pd.read_csv(f'{DATA}/fig3c_tau.csv')
BM  = pd.read_csv(f'{DATA}/fig3d_benchmark.csv')

fig = plt.figure(figsize=(15.5, 12.0))
gs = fig.add_gridspec(2, 2, hspace=0.58, wspace=0.30, left=0.07, right=0.975)

# (a) tier enrichment
ax = fig.add_subplot(gs[0, 0])
EN = EN.iloc[::-1].reset_index(drop=True)
yp = np.arange(len(EN)); off = 0.18
ax.axvline(1, color='black', lw=1.1)
ax.scatter(EN['sink_OR'], yp+off, s=150, color=BLUE, zorder=3,
           edgecolors='white', linewidths=1.3, label='sink tier (15th pct)')
ax.scatter(EN['source_OR'], yp-off, s=150, marker='s', color=ORANGE, zorder=3,
           edgecolors='white', linewidths=1.3, label='source tier (85th pct)')
ax.set_yticks(yp)
ax.set_yticklabels([s.replace(' (', '\n(') for s in EN['variant_set']], fontsize=FS['tick']-1)
ax.set_xscale('log'); ax.set_xlim(0.15, 6)
ax.set_xticks([0.2, 0.5, 1, 2, 5]); ax.set_xticklabels(['0.2', '0.5', '1', '2', '5'])
ax.xaxis.set_minor_locator(mpl.ticker.NullLocator())
ax.set_xlabel('Odds ratio (NBD-dimerized consensus)')
ax.legend(frameon=False, loc='lower right')
ax.spines['left'].set_visible(False); ax.tick_params(axis='y', length=0)
panel(ax, 'a')

# (b) netTE vs AlphaMissense
ax = fig.add_subplot(gs[0, 1])
ax.scatter(SC['netTE_z'], SC['AM_balance'], s=14, alpha=0.30, color=BLUE, edgecolors='none')
z_ = np.polyfit(SC['netTE_z'], SC['AM_balance'], 1)
xx = np.linspace(SC['netTE_z'].min(), SC['netTE_z'].max(), 50)
ax.plot(xx, np.polyval(z_, xx), color=RED, lw=2.4)
ax.set_xlabel('netTE (z, NBD-dimerized consensus)')
ax.set_ylabel('AlphaMissense balance\n(pathogenic $-$ benign substitutions)')
ax.text(0.03, 0.05, '$\\rho$ = $-$0.560\nn = %d' % len(SC), transform=ax.transAxes,
        fontsize=FS['annot'], va='bottom')
panel(ax, 'b')

# (c) tau dependence of the source/sink assignment
ax = fig.add_subplot(gs[1, 0])
for pid, mk, col, lab in [('9dw9', 'o', BLUE, '9DW9 (dimerized)'),
                          ('5uak', 's', ORANGE, '5UAK (separated)')]:
    s_ = TAU[TAU['structure'] == pid].sort_values('tau')
    ax.plot(s_['tau'], s_['sink_OR'], marker=mk, color=col, lw=2.2, ms=9, label=f'{lab}, sink')
    ax.plot(s_['tau'], s_['source_OR'], marker=mk, color=col, lw=2.2, ms=9, ls='--',
            markerfacecolor='white', label=f'{lab}, source')
ax.axhline(1, color='black', lw=1.1)
ax.set_xscale('log'); ax.set_yscale('log')
ax.set_yticks([0.3, 0.5, 1, 2, 3]); ax.set_yticklabels(['0.3', '0.5', '1', '2', '3'])
ax.yaxis.set_minor_locator(mpl.ticker.NullLocator())
ax.xaxis.set_minor_locator(mpl.ticker.NullLocator())
ax.set_xlabel('transfer entropy time delay $\\tau$')
ax.set_ylabel('Odds ratio, CF-causing variants')
ax.legend(frameon=False, fontsize=FS['legend']-3, ncol=2, loc='upper center',
          bbox_to_anchor=(0.5, 1.26))
panel(ax, 'c')

# (d) discrimination benchmark
ax = fig.add_subplot(gs[1, 1])
colmap = {'residue': '#7fa8cc', 'substitution': BLUE}
xs = np.arange(len(BM))
for i, (_, r) in enumerate(BM.iterrows()):
    c = GREEN if 'ceiling' in r['label'] else colmap[r['level']]
    ax.bar(i, r['AUC'], width=0.68, color=c, edgecolor='white', zorder=2)
    ax.text(i, r['AUC']+0.008, '%.3f' % r['AUC'], ha='center', fontsize=FS['annot']-1)
ax.axvline(2.5, color=GREY, lw=1.2, ls='--')
ax.set_xticks(xs)
ax.set_xticklabels([l.replace('netTE + burial\n+ chemistry', 'netTE\n+burial\n+chem')
                     .replace('per-residue\noracle ceiling', 'oracle\nceiling') for l in BM['label']],
                   fontsize=FS['tick']-2)
ax.set_ylim(0.5, 0.98); ax.set_ylabel('AUC')
ax.text(1.0, 0.955, 'residue level\n(CF-causing vs not)', ha='center',
        fontsize=FS['annot']-1, color='0.35')
ax.text(4.0, 0.955, 'substitution level\n(cross-validated)', ha='center',
        fontsize=FS['annot']-1, color='0.35')
panel(ax, 'd')

fig.savefig(f'{OUT}/figure3_scope_severity.png', dpi=300, bbox_inches='tight')
fig.savefig(f'{OUT}/figure3_scope_severity.pdf', bbox_inches='tight')
plt.close(fig)

# ============================ FIGURE 3 =======================================
d  = pd.read_csv(f'{DATA}/fig4b_scatter.csv')
gr = pd.read_csv(f'{DATA}/fig4c_grid_rho.csv', index_col=0)
gp = pd.read_csv(f'{DATA}/fig4c_grid_p.csv',  index_col=0)
F  = pd.read_csv(f'{DATA}/fig4d_forest.csv')
DZ = pd.read_csv(f'{DATA}/deltaz_consensus.csv', index_col=0)['delta_z']; DZ.index = DZ.index.astype(int)
LAB = pd.read_csv(f'{DATA}/key_variant_labels.csv')
DOM = [(1,390,'TMD1'),(390,650,'NBD1'),(650,860,'R'),(860,1170,'TMD2'),(1170,1480,'NBD2')]

fig = plt.figure(figsize=(15.5, 15.0))
gs = fig.add_gridspec(3, 2, height_ratios=[0.95, 1.0, 0.92], width_ratios=[1.0, 1.45],
                      hspace=0.42, wspace=0.30, left=0.06, right=0.965)

# (a) Delta z along the sequence
ax = fig.add_subplot(gs[0, :])
# reserved label bands: positives above the trace maximum, negatives below the minimum,
# so no annotation can land on the domain row or the x-axis
TIER = 0.78
hi_base, lo_base = DZ.max() + 0.75, DZ.min() - 0.75
y_hi = hi_base + 2 * TIER          # top of the positive label band
y_lo = lo_base - 2 * TIER          # bottom of the negative label band
ytop = y_hi + 1.75                 # domain-name row, clear of every label
ax.set_xlim(1, 1480); ax.set_ylim(y_lo - 0.85, ytop + 0.45)

for x0, x1, nm in DOM:
    if nm in ('NBD1', 'NBD2'):
        ax.axvspan(x0, x1, color='#f4f4f4', zorder=0)
    ax.axvline(x0, color='0.86', lw=0.9, zorder=0)
    ax.text((x0+x1)/2, ytop, nm, ha='center', va='top', fontsize=FS['annot'], color='0.35')
ax.axhline(0, color='black', lw=1.0, zorder=1)
ax.fill_between(DZ.index, 0, DZ.values, where=DZ.values > 0, color=BLUE, alpha=0.30, interpolate=True, zorder=2)
ax.fill_between(DZ.index, 0, DZ.values, where=DZ.values < 0, color=ORANGE, alpha=0.30, interpolate=True, zorder=2)
ax.plot(DZ.index, DZ.values, lw=0.9, color='0.25', zorder=3)

_lab = LAB[LAB['residue'].astype(int).isin(DZ.index)].copy()
_lab['residue'] = _lab['residue'].astype(int)
_lab['dz'] = _lab['residue'].map(DZ)
for sign in (1, -1):
    sub = _lab[np.sign(_lab['dz']) == sign].sort_values('residue').reset_index(drop=True)
    for k, r in sub.iterrows():
        rr, y = int(r['residue']), r['dz']
        ax.scatter([rr], [y], s=62, color=RED, zorder=5, edgecolors='white', linewidths=1.0)
        tier = k % 3
        ylab = (hi_base + tier * TIER) if sign > 0 else (lo_base - tier * TIER)
        ax.annotate(r['variant'], xy=(rr, y), xytext=(rr, ylab), textcoords='data',
                    ha='center', va='bottom' if sign > 0 else 'top',
                    fontsize=FS['annot'] - 2, color=RED,
                    arrowprops=dict(arrowstyle='-', color=RED, lw=0.8, alpha=0.7,
                                    shrinkA=1.5, shrinkB=3))
ax.set_xlabel('CFTR residue number')
ax.set_ylabel('$\\Delta$z')
panel(ax, 'a')

# (b) scatter
ax = fig.add_subplot(gs[1, 0])
ax.axhline(0, color=GREY, lw=1.0, ls=':', zorder=1)
ax.scatter(d['dz'], d['corrector_gain_CB'], s=20, alpha=0.35, color=BLUE, edgecolors='none', zorder=2)
q = pd.qcut(d['dz'], 5, labels=False); rb = np.random.default_rng(3)
cx, cy, lo, hi = [], [], [], []
for k in range(5):
    sb = d[q == k]
    cx.append(sb['dz'].median()); cy.append(sb['corrector_gain_CB'].median())
    bs = [np.median(rb.choice(sb['corrector_gain_CB'].values, len(sb))) for _ in range(2000)]
    lo.append(np.percentile(bs, 2.5)); hi.append(np.percentile(bs, 97.5))
ax.errorbar(cx, cy, yerr=[np.array(cy)-np.array(lo), np.array(hi)-np.array(cy)], fmt='o-',
            color=RED, ms=9, lw=2.4, capsize=4, zorder=4, markeredgecolor='white',
            markeredgewidth=1.2, label='quintile median (95% CI)')
ax.set_ylim(-160, 300)
ax.set_xlabel('$\\Delta$z')
ax.set_ylabel('Corrector gain, C/B ratio\n(ELX/TEZ $-$ DMSO, % WT)')
ax.legend(frameon=False, loc='upper left')
ax.text(0.97, 0.05, '$\\rho$ = %.3f\nP = 1.5$\\times$10$^{-11}$\nn = %d' % (F.loc[0,'rho'], F.loc[0,'n']),
        transform=ax.transAxes, ha='right', va='bottom', fontsize=FS['annot'])
panel(ax, 'b')

# (c) predictor x measurement grid
ax = fig.add_subplot(gs[1, 1])
Gv, Pv = gr.values, gp.values
im = ax.imshow(Gv, cmap='RdBu_r', vmin=-0.45, vmax=0.45, aspect='auto')
SHORT = {'Corrector gain (C/B ratio)':'Corrector\ngain (C/B)','Corrector gain (band C)':'Corrector\ngain (band C)',
         'Function gain (FSK)':'Function\ngain (FSK)','Combined ELX/TEZ/IVA response':'Combined\nELX/TEZ/IVA',
         'Potentiator gain (IVA)':'Potentiator\ngain (IVA)'}
ax.set_xticks(range(Gv.shape[1]))
ax.set_xticklabels([SHORT.get(c, c) for c in gr.columns], fontsize=FS['tick']-1,
                   rotation=22, ha='right', rotation_mode='anchor')
ax.set_yticks(range(Gv.shape[0]))
ax.set_yticklabels(['$\\Delta$z' if i.startswith('\u0394z') else ('AlphaMissense' if 'Alpha' in i else 'Burial')
                    for i in gr.index], fontsize=FS['tick'])
for i in range(Gv.shape[0]):
    for jj in range(Gv.shape[1]):
        stt = '***' if Pv[i,jj] < 1e-4 else ('**' if Pv[i,jj] < 1e-2 else ('*' if Pv[i,jj] < 0.05 else ''))
        ax.text(jj, i, '%+.2f\n%s' % (Gv[i,jj], stt), ha='center', va='center',
                fontsize=FS['heat_cell'], color='white' if abs(Gv[i,jj]) > 0.25 else 'black')
ax.add_patch(Rectangle((-0.5,-0.5), 1, Gv.shape[0], fill=False, edgecolor=GREEN, lw=3.0, zorder=5))
ax.add_patch(Rectangle((Gv.shape[1]-1.5,-0.5), 1, Gv.shape[0], fill=False, edgecolor=ORANGE, lw=2.6, ls='--', zorder=5))
cb = fig.colorbar(im, ax=ax, fraction=0.030, pad=0.02)
cb.set_label('Spearman $\\rho$', fontsize=FS['axis_label']); cb.ax.tick_params(labelsize=FS['tick'])
panel(ax, 'c')

# (d) forest
ax = fig.add_subplot(gs[2, :])
yp = np.arange(len(F))[::-1]; off = 0.17
ax.axvline(0, color='black', lw=1.1, zorder=1)
ax.errorbar(F['rho'], yp+off, xerr=[F['rho']-F['lo'], F['hi']-F['rho']], fmt='o', color=BLUE,
            ms=11, lw=2.6, capsize=5, zorder=3, markeredgecolor='white', markeredgewidth=1.2,
            label='$\\Delta$z')
ax.scatter(F['rho_am'], yp-off, marker='s', s=90, color=ORANGE, zorder=3,
           edgecolors='white', linewidths=1.2, label='AlphaMissense')
for i, (_, r) in enumerate(F.iterrows()):
    ax.text(0.70, yp[i]+off, '%+.3f  [%.2f, %.2f]' % (r['rho'], r['lo'], r['hi']),
            fontsize=FS['annot'], va='center', color=BLUE)
    ax.text(0.70, yp[i]-off, '%+.3f' % r['rho_am'], fontsize=FS['annot'], va='center', color=ORANGE)
ax.set_yticks(yp)
ax.set_yticklabels(['%s\n%s  (n = %d)' % (r['dataset'], r['assay'], r['n']) for _, r in F.iterrows()],
                   fontsize=FS['ytick_forest'])
ax.set_ylim(-0.62, len(F)-0.18); ax.set_xlim(-0.33, 1.00)
ax.set_xticks(np.arange(-0.3, 0.75, 0.1))
ax.set_xlabel('Spearman $\\rho$ with corrector-induced rescue  (95% CI, residue-level bootstrap)')
ax.legend(frameon=False, loc='upper center', ncol=2, bbox_to_anchor=(0.45, 1.10))
ax.spines['left'].set_visible(False); ax.tick_params(axis='y', length=0)
panel(ax, 'd')

fig.savefig(f'{OUT}/figure4_deltaz_corrector_rescue.png', dpi=300, bbox_inches='tight')
fig.savefig(f'{OUT}/figure4_deltaz_corrector_rescue.pdf', bbox_inches='tight')
plt.close(fig)

# ============================ FIGURE 4 (was 3) ===============================
NC = pd.read_csv(f'{DATA}/fig5b_within_state_controls.csv')
NULL = pd.read_csv(f'{DATA}/fig5c_permutation_null.csv')
OBS  = pd.read_csv(f'{DATA}/fig5c_permutation_observed.csv')
NS   = pd.read_csv(f'{DATA}/fig5d_truncating.csv')
CP   = pd.read_csv(f'{DATA}/fig5a_drugclass.csv')

fig = plt.figure(figsize=(15.5, 12.0))
gs = fig.add_gridspec(2, 2, hspace=0.52, wspace=0.24)
cols = {'bihler': BLUE, 'tedman': '#5aa2d8', 'ajrccm': '#a9cce6'}

ax = fig.add_subplot(gs[0, 1])
rows_ = NC[NC['kind'] != 'permutation null mean'].reset_index(drop=True)
labs = ['9MXL $-$ 6MSM\nwithin dimerized', '9DW9 $-$ 8V81\nwithin dimerized',
        '5UAK $-$ 9DW5\nwithin separated', 'half A $-$ half B\nwithin dimerized',
        'dimerized $-$ separated\nREAL $\\Delta$z']
xs = np.arange(len(rows_)); w = 0.26
ax.axhspan(-0.17, 0.17, color='#eeeeee', zorder=0)
ax.axhline(0, color='black', lw=1.0, zorder=1)
for k, (c, lab) in enumerate([('bihler', 'Bihler'), ('tedman', 'Tedman'), ('ajrccm', 'Castellani')]):
    ax.bar(xs+(k-1)*w, rows_[c], width=w, color=cols[c], label=lab, zorder=2,
           edgecolor='white', linewidth=0.8)
ax.axvline(3.5, color=GREY, lw=1.2, ls='--')
ax.set_xticks(xs); ax.set_xticklabels(labs, fontsize=FS['tick']-2, rotation=25,
                                      ha='right', rotation_mode='anchor')
ax.set_ylabel('Spearman $\\rho$ with corrector rescue')
ax.legend(frameon=False, loc='upper left', ncol=3, columnspacing=0.9)
ax.set_ylim(-0.24, 0.50)
panel(ax, 'b')

ax = fig.add_subplot(gs[1, 0])
pc = {'Bihler': BLUE, 'Tedman': '#5aa2d8', 'AJRCCM': '#a9cce6'}
for c in NULL.columns:
    ax.hist(NULL[c].values, bins=40, histtype='stepfilled', alpha=0.45,
            color=pc.get(c, BLUE), edgecolor=pc.get(c, BLUE), lw=1.2, label=f'{c} null')
ymax = ax.get_ylim()[1]; ax.set_ylim(0, ymax*1.38)
heights = {'Castellani': 0.40, 'Bihler': 0.64, 'Tedman': 0.40}
for _, r in OBS.iterrows():
    h = heights.get(r['dataset'], 0.5); col = pc.get(r['dataset'], BLUE)
    ax.annotate('', xy=(r['observed_rho'], ymax*0.04), xytext=(r['observed_rho'], ymax*h),
                arrowprops=dict(arrowstyle='-|>', color=col, lw=2.4))
    ax.text(r['observed_rho'], ymax*(h+0.03), '%s\n$\\rho$ = %.3f' % (r['dataset'], r['observed_rho']),
            fontsize=FS['annot']-1, ha='center', color=col)
ax.set_xlim(-0.25, 0.52)
ax.set_xlabel('Spearman $\\rho$ under residue-label permutation (1,000$\\times$)')
ax.set_ylabel('Permutations')
ax.legend(frameon=False, loc='upper left')
ax.text(0.99, 0.99, 'all three: empirical P < 0.001\n%.1f$-$%.1f null SD from mean'
        % (OBS['sd_from_null'].min(), OBS['sd_from_null'].max()),
        transform=ax.transAxes, ha='right', va='top', fontsize=FS['annot'])
panel(ax, 'c')

ax = fig.add_subplot(gs[1, 1])
ax.axhline(0, color='black', lw=1.0)
for i, (_, r) in enumerate(NS.iterrows()):
    ax.bar(i, r['rho'], width=0.52, color=[BLUE, GREY][i], edgecolor='white', zorder=2)
    ax.errorbar(i, r['rho'], yerr=[[r['rho']-r['lo']], [r['hi']-r['rho']]], fmt='none',
                ecolor='black', capsize=5, lw=1.8, zorder=3)
    ax.text(i, r['hi']+0.05, '$\\rho$ = %+.3f\nP = %.2g\nn = %d' % (r['rho'], r['p'], r['n']),
            ha='center', fontsize=FS['annot'])
ax.set_xticks([0, 1])
ax.set_xticklabels(['Missense\n(signal expected)', 'Truncating\n(no folded protein)'])
ax.set_ylabel('$\\rho$, $\\Delta$z vs corrector fold-change')
ax.set_ylim(-0.88, 1.02)
panel(ax, 'd')

ax = fig.add_subplot(gs[0, 0])
ax.axhline(0, color='black', lw=1.0)
outs = ['Corrector gain', 'Potentiator gain']; w = 0.32
for k, (pl, col) in enumerate([('\u0394z', BLUE), ('AlphaMissense', ORANGE)]):
    key = CP['predictor'].unique()[k]
    sub = CP[CP['predictor'] == key].set_index('outcome').loc[outs]
    lab = '$\\Delta$z' if k == 0 else 'AlphaMissense'
    ax.bar(np.arange(2)+(k-0.5)*w, sub['rho'], width=w, color=col, label=lab,
           edgecolor='white', linewidth=0.8, zorder=2)
    ax.errorbar(np.arange(2)+(k-0.5)*w, sub['rho'],
                yerr=[sub['rho']-sub['lo'], sub['hi']-sub['rho']], fmt='none',
                ecolor='black', capsize=5, lw=1.6, zorder=3)
ax.set_xticks([0, 1])
ax.set_xticklabels(['Corrector rescue\n(ELX/TEZ, trafficking)', 'Potentiator rescue\n(IVA, gating)'])
ax.set_ylabel('Spearman $\\rho$ (95% CI)')
ax.set_ylim(-0.42, 0.50)
ax.legend(frameon=False, loc='upper right')
panel(ax, 'a')

fig.savefig(f'{OUT}/figure5_specificity_controls.png', dpi=300, bbox_inches='tight')
fig.savefig(f'{OUT}/figure5_specificity_controls.pdf', bbox_inches='tight')
plt.close(fig)
print('wrote figures to', os.path.normpath(OUT))
