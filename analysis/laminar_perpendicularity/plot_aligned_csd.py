"""Sound-aligned CSD for selected probes.
Aligns to every element onset (5 per sequence) = stimulusOnset_ms + per-trial audio latency
(Walt: measured from the recorded audio, session_audio_latency.mat; Troy: fixed 200 ms placeholder - no audio recorded).
Same preprocessing as laminar_perpendicularity.py (bad-contact interpolation, gain equalisation).
Usage: python3 plot_aligned_csd.py <results_dir> session:probe [session:probe ...]"""
import sys, os, json, numpy as np, h5py, scipy.io as sio, openpyxl, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt, seaborn as sns
import laminar_perpendicularity as lp
sns.set_theme(style='ticks', context='paper')
REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
DATA = '/Volumes/Mnemosyne/Data/2025_macaque_aglt' if os.path.isdir('/Volumes/Mnemosyne/Data/2025_macaque_aglt') else os.path.expanduser('~/mnt/2025_macaque_aglt')
ELEM = np.array([0, 563, 1126, 1689, 2252])
PRE, POST = 100, 400

def latency_map():
    """session -> per-trial audio latency (ms). Matches rows of session_audio_latency.mat to sessions by trial count
    (the .mat was built from a log version missing one session, so row order drifts by one after row 61)."""
    wb = openpyxl.load_workbook(f'{REPO}/data-extraction/doc/log_backup/kikuchi_recording_log_20240812.xlsx', read_only=True, data_only=True)
    ses = [r[5] for r in list(wb['agl_t'].iter_rows(values_only=True))[1:] if r[5]][0::2]
    lat = sio.loadmat(f'{REPO}/data-extraction/session_audio_latency.mat')['session_audio_latency']
    out, j = {}, 0
    for i in range(lat.shape[0]):
        x = np.asarray(lat[i, 0], float).ravel()
        while j < len(ses):
            with h5py.File(f'{DATA}/{ses[j]}.mat', 'r') as f: m = lp.read_event_table(f)['trial_n'].size
            j += 1
            if m == x.size: out[ses[j - 1]] = x; break
    return out

def main():
    res = sys.argv[1]; targets = [a.split(':') for a in sys.argv[2:]]
    lm = latency_map()
    fig, ax = plt.subplots(len(targets), 3, figsize=(12, 3.1 * len(targets)), gridspec_kw=dict(width_ratios=[2.2, 1.4, 0.9]), squeeze=False)
    for i, (s, p) in enumerate(targets):
        p = int(p); meta = [r for r in json.load(open(f'{res}/probe_data/{s}.json')) if r['probe'] == p][0]
        with h5py.File(f'{DATA}/{s}.mat', 'r') as f:
            ev = lp.read_event_table(f); lfp = f['lfp'][()].astype(np.float32)
        ch0 = (p - 1) * 16; x = lfp[:, ch0:ch0 + 16]
        x = lp.patch(x, lp.detect_bad(x)); x, _ = lp.equalise_gain(x)
        on = ev['stimulusOnset_ms']
        if s.startswith('walt') and s in lm:
            lat = lm[s].copy(); lat[(lat < 0) | (lat > 300)] = np.nan; src = 'measured audio latency'
        else:
            lat = np.full(on.shape, 200.); src = 'fixed 200 ms (no audio)'
        t0 = on + lat; t0 = t0[np.isfinite(t0)]
        el = (t0[:, None] + ELEM[None, :]).ravel()
        ep = lp.epochs(x, el, PRE, POST).astype(float); ep -= ep[:, :PRE].mean(1, keepdims=True)
        erp = ep.mean(0); c = lp.csd(erp, meta['spacing_um'])
        c_odd, c_even = lp.csd(ep[0::2].mean(0), meta['spacing_um']), lp.csd(ep[1::2].mean(0), meta['spacing_um'])
        w = slice(PRE, PRE + 150); rsh = np.corrcoef(c_odd[w].ravel(), c_even[w].ravel())[0, 1]
        tt = np.arange(-PRE, POST)
        # CSD
        a = ax[i, 0]; v = np.nanpercentile(np.abs(c[PRE:PRE + 300]), 99)
        im = a.imshow(c.T, aspect='auto', cmap='RdBu', vmin=-v, vmax=v, extent=[-PRE, POST, 16.5, .5])
        a.axvline(0, color='k', lw=.6); a.axhline(meta['crossover_ch'], color='k', ls='--', lw=.8)
        a.set(ylabel='contact', title=f"{s}  P{p} ({meta['area']})  n={len(ep)} elements\nCSD (red = source, blue = sink); onset = {src}; split-half r={rsh:.2f}")
        if i == len(targets) - 1: a.set_xlabel('ms from sound element onset')
        # ERP traces
        a = ax[i, 1]; sc = np.nanmax(np.abs(erp[PRE:])) * 0.6 + 1e-9
        for k in range(16): a.plot(tt, -k + erp[:, k] / sc, color='k', lw=.7)
        a.axvline(0, color='k', lw=.6); a.set(yticks=-np.arange(0, 16, 3), yticklabels=np.arange(1, 17, 3), title='evoked LFP (gain-equalised)')
        if i == len(targets) - 1: a.set_xlabel('ms')
        # spectral profile
        d = np.load(f'{res}/probe_data/{s}_p{p}.npz'); chn = np.arange(1, 17)
        a = ax[i, 2]; a.plot(d['ab'], chn, color='#ff7f0e', label='α-β'); a.plot(d['g'], chn, color='#9467bd', label='γ')
        a.axhline(meta['crossover_ch'], color='k', ls='--', lw=.8); a.set_ylim(16.5, .5); a.set(title=f"G={meta['G']:.2f}", xlabel='rel. power')
        a.legend(fontsize=7, frameon=False)
    fig.tight_layout(); out = f'{res}/aligned_csd_perpendicular.png'; fig.savefig(out, dpi=120); print(out)

if __name__ == '__main__':
    main()
