"""One QC panel row per probe: relative power profiles, gamma-band correlation matrix, evoked CSD.
Usage: python3 plot_probe_qc.py <out_dir> [session ...]  -> qc_figures/<session>.png"""
import sys, os, json, glob, numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt, seaborn as sns
sns.set_theme(style='ticks', context='paper')
out = sys.argv[1]; os.makedirs(os.path.join(out, 'qc_figures'), exist_ok=True)
sess = sys.argv[2:] or sorted(os.path.basename(f)[:-5] for f in glob.glob(os.path.join(out, 'probe_data', '*.json')))
AREA_COL = {'auditory': '#1f77b4', 'frontal': '#d62728', 'hpc': '#2ca02c'}
for s in sess:
    rows = json.load(open(os.path.join(out, 'probe_data', s + '.json')))
    fig, ax = plt.subplots(2, 4, figsize=(13, 6.5), gridspec_kw=dict(width_ratios=[1.3, 1, 1, 1]))
    for i, r in enumerate(rows):
        d = np.load(os.path.join(out, 'probe_data', f"{s}_p{r['probe']}.npz"))
        ch = np.arange(1, 17)
        a = ax[i, 0]; im = a.imshow(d['rel'][:, d['f'] <= 150], aspect='auto', cmap='viridis',
                                    extent=[0, 150, 16.5, 0.5]); a.set(xlabel='Hz', ylabel='contact', title=f"P{r['probe']} {r['area']}: rel. power")
        a = ax[i, 1]; a.plot(d['ab'], ch, color='#ff7f0e', label='α-β 10–19'); a.plot(d['g'], ch, color='#9467bd', label='γ 75–145')
        if 'crossover_ch' in r and np.isfinite(r['crossover_ch']): a.axhline(r['crossover_ch'], color='k', ls='--', lw=.8)
        a.invert_yaxis(); a.set(xlabel='relative power', title=f"G={r['G']:.2f} p={r['p_motif']:.3f}"); a.legend(fontsize=7, frameon=False)
        a = ax[i, 2]; a.imshow(d['gammaC'], vmin=0, vmax=1, cmap='mako', extent=[.5, 16.5, 16.5, .5])
        a.axhline(r['block_split_after'] + .5, color='w', lw=.8); a.set(title=f"γ corr; block={r['block_index']:.2f}")
        a = ax[i, 3]; c = d['csd'][100:300]; v = np.nanpercentile(np.abs(c), 98) if np.isfinite(c).any() else 1
        a.imshow(c.T, aspect='auto', cmap='RdBu', vmin=-v, vmax=v, extent=[0, 200, 16.5, .5])
        a.set(xlabel='ms from element onset', title=f"CSD snr={r['evoked_snr']:.1f} r½={r['csd_splithalf_r']:.2f}")
        for a in ax[i]: [sp.set_color(AREA_COL.get(r['area'], 'gray')) for sp in a.spines.values()]
    fig.suptitle(f"{s}   P1 ({rows[0]['area']}): {rows[0]['verdict']} [{rows[0].get('confidence','')}]   |   P2 ({rows[1]['area']}): {rows[1]['verdict']} [{rows[1].get('confidence','')}]", fontsize=10)
    fig.tight_layout(); fig.savefig(os.path.join(out, 'qc_figures', s + '.png'), dpi=110); plt.close(fig)
