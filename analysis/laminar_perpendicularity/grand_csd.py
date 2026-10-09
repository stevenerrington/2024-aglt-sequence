"""Grand-average sound-evoked CSD per area, depth-aligned to each probe's dominant evoked current sink.

Step 1 (compute, cached per probe in <results>/aligned_csd/): element-onset-aligned CSD (5 elements per sequence),
  onset = stimulusOnset_ms + audio latency (Walt: measured per trial, session median if missing / missing session
  -> monkey median; Troy: fixed 200 ms placeholder, no audio recorded). Same preprocessing as laminar_perpendicularity.py.
Step 2 (plot): orientation set so superficial is up (gamma-peak end: Troy low index, Walt high index); depth in um;
  sink = most negative CSD in 0-200 ms (mean over 20-ms window around its time), required z < -3 vs baseline
  (-100-0 ms); each probe normalised by its peak |CSD| (0-300 ms), interpolated onto a 50-um grid relative to the
  sink, averaged per area x monkey. CSD sign: negative = sink (blue), positive = source (red).
Usage: python3 grand_csd.py <results_dir> compute|plot"""
import sys, os, json, glob, numpy as np, h5py, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt, seaborn as sns, pandas as pd
import laminar_perpendicularity as lp
from plot_aligned_csd import latency_map, DATA, ELEM
PRE, POST = 100, 400
res, mode = sys.argv[1], sys.argv[2]
outd = f'{res}/aligned_csd'; os.makedirs(outd, exist_ok=True)
S = pd.read_csv(f'{res}/perpendicularity_summary.csv')
S = S[S.area_group.isin(['auditory', 'frontal'])]

if mode == 'compute':
    lm = None
    for s in sorted(S.session.unique()):
        if os.path.exists(f'{outd}/{s}.npz'): continue
        if lm is None: lm = latency_map(); walt_med = np.nanmedian(np.concatenate([v for k, v in lm.items() if k.startswith('walt')]))
        with h5py.File(f'{DATA}/{s}.mat', 'r') as f:
            ev = lp.read_event_table(f); lfp = f['lfp'][()].astype(np.float32)
        on = ev['stimulusOnset_ms']
        if s.startswith('walt'):
            lat = lm[s].copy() if s in lm else np.full(on.shape, walt_med)
            bad = ~np.isfinite(lat) | (lat < 0) | (lat > 300); lat[bad] = np.nanmedian(lat[~bad]) if (~bad).any() else walt_med
        else:
            lat = np.full(on.shape, 200.)
        t0 = on + lat; t0 = t0[np.isfinite(t0)]; el = (t0[:, None] + ELEM[None, :]).ravel()
        out = {}
        for p in (1, 2):
            row = S[(S.session == s) & (S.probe == p)]
            if row.empty: continue
            meta = [r for r in json.load(open(f'{res}/probe_data/{s}.json')) if r['probe'] == p][0]
            x = lfp[:, (p - 1) * 16:p * 16]; x = lp.patch(x, lp.detect_bad(x)); x, _ = lp.equalise_gain(x)
            ep = lp.epochs(x, el, PRE, POST).astype(float); ep -= ep[:, :PRE].mean(1, keepdims=True)
            c = lp.csd(ep.mean(0), meta['spacing_um'])
            # baseline noise of CSD from trial-shuffled halves
            out[f'p{p}'] = c; out[f'p{p}_base_sd'] = c[:PRE].std(0).mean(); out[f'p{p}_spacing'] = meta['spacing_um']
        np.savez_compressed(f'{outd}/{s}.npz', **out); print(s, flush=True)

if mode == 'plot':
    sns.set_theme(style='ticks', context='paper')
    grid = np.arange(-1500, 1501, 50); groups = {}; rows = []
    for _, r in S.iterrows():
        z = np.load(f'{outd}/{r.session}.npz'); c = z[f'p{r.probe}']; sp = float(z[f'p{r.probe}_spacing'])
        if not np.isfinite(c).all(): continue
        if r.monkey == "walt": c = c[:, ::-1]            # superficial up
        base = c[:PRE]; w = c[PRE:PRE + 200]
        tmin, kmin = np.unravel_index(np.argmin(w), w.shape)
        sink_val = w[max(0, tmin - 10):tmin + 10, kmin].mean()
        zsink = (sink_val - base[:, kmin].mean()) / (base.std() + 1e-12)
        rows.append(dict(session=r.session, probe=r.probe, monkey=r.monkey, area=r.area_group, verdict=r.verdict,
                         sink_contact=kmin + 1, sink_ms=tmin, sink_z=zsink))
        if zsink > -3: continue
        cn = c / (np.abs(c[PRE:PRE + 300]).max() + 1e-12)
        depth = (np.arange(16) - kmin) * sp
        ci = np.stack([np.interp(grid, depth, cn[t], left=np.nan, right=np.nan) for t in range(cn.shape[0])])
        for key in [(r.area_group, r.monkey, 'all'), (r.area_group, r.monkey, 'perp/probable') if r.verdict in ('perpendicular', 'probable') else None]:
            if key: groups.setdefault(key, []).append(ci)
    T = pd.DataFrame(rows); T.to_csv(f'{res}/sink_alignment.csv', index=False)
    print(T.groupby(['area', 'monkey']).apply(lambda d: f"{(d.sink_z < -3).sum()}/{len(d)} with sink z<-3; median sink t={d[d.sink_z<-3].sink_ms.median()} ms").to_string())
    keys = [(a, m, v) for a in ('auditory', 'frontal') for m in ('troy', 'walt') for v in ('all', 'perp/probable')]
    fig, ax = plt.subplots(2, 4, figsize=(16, 8), sharey=True)
    tt = np.arange(-PRE, POST)
    for i, (a, m, v) in enumerate(keys):
        A = ax[i // 4, i % 4]; L = groups.get((a, m, v), [])
        if not L: A.set_title(f'{a} | {m} | {v}: none'); continue
        G = np.nanmean(np.stack(L), 0); n = np.sum(np.isfinite(np.stack(L)[:, PRE + 50, :]), 0)
        G[:, n < max(3, len(L) // 3)] = np.nan
        vmax = np.nanpercentile(np.abs(G[PRE:PRE + 300]), 99)
        A.imshow(G.T, aspect='auto', cmap='RdBu_r', vmin=-vmax, vmax=vmax, extent=[-PRE, POST, grid[-1] + 25, grid[0] - 25])
        A.axhline(0, color='k', lw=.6, ls='--'); A.axvline(0, color='k', lw=.6)
        A.set(title=f'{a} | {m} | {v}  (n={len(L)})', xlim=(-PRE, 300))
        if i % 4 == 0: A.set_ylabel('µm relative to evoked sink\n(− = superficial)')
        if i // 4 == 1: A.set_xlabel('ms from sound element onset')
    fig.suptitle('Grand-average sound-evoked CSD, aligned to each probe\'s evoked sink (blue = sink, red = source). '
                 'Troy onsets use fixed 200 ms latency (no audio); Walt onsets use measured audio latency.', fontsize=10)
    fig.tight_layout(); fig.savefig(f'{res}/grand_csd_sink_aligned.png', dpi=120); print('saved')
