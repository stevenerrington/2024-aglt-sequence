"""
Per-probe test of whether each AGL-t penetration ran perpendicular to the cortical layers.

Three independent signatures, computed per 16-contact probe (contacts 1-16 = probe 1, 17-32 = probe 2):
  1. Spectrolaminar motif (Mendoza-Halliday et al., 2024, Nat Neurosci): relative alpha-beta (10-19 Hz)
     power increases with depth while relative gamma (75-145 Hz, line-noise harmonics removed) decreases,
     crossing near L4. Fit within the best contiguous window of >= 8 contacts; goodness G = R2_ab * R2_gamma
     when slopes have opposite signs (else 0). Significance from shuffling contact order (1000 perms).
  2. Evoked-LFP polarity reversal: the sound-evoked LFP (0-150 ms) inverts across depth in a perpendicular
     penetration; a tangential one gives near-identical waveforms on all contacts.
  0. Preprocessing: data-driven bad-contact interpolation; contact gain equalisation (equal 1-8 Hz RMS) because TDT/Plexon contacts differ in gain/impedance.
  2b. Gamma-band (30-145 Hz) inter-contact correlation block structure; counts only if its split lies
     within 2 contacts of the spectral crossover.
  3. Evoked CSD reliability: split-half (odd vs even trials) correlation of the CSD map (0-150 ms).
     Perpendicular penetrations produce a stable sink/source pattern; tangential ones give noise.

Self-contained: reads MATLAB v7.3 session files (lfp at 1 kHz, event_table) directly.
Usage: python3 laminar_perpendicularity.py <mat_dir> <log_xlsx> <out_dir> [session ...]
"""
import sys, os, json, numpy as np, h5py, openpyxl
from scipy.signal import welch
from scipy.stats import linregress

FS = 1000
AB_BAND = (10, 19)
G_BANDS = [(75, 95), (105, 145)]          # skip 100 Hz harmonic of UK mains
MIN_WIN = 8
N_PERM = 1000
rng = np.random.default_rng(0)

# --------------------------------------------------------------------------- I/O
def read_event_table(f):
    m = f['#subsystem#/MCOS'][()].ravel()
    objs = [f[r] for r in m if isinstance(f[r], h5py.Dataset) and f[r].dtype == object]
    for d, v in zip(objs[:-1], objs[1:]):
        try:
            names = [''.join(map(chr, f[x][()].ravel())) for x in v[()].ravel()]
        except Exception:
            continue
        if len(names) == d.size and all(n.isidentifier() for n in names):
            return {n: (f[x][()].ravel() if f[x].dtype.kind in 'fiu' else None)
                    for n, x in zip(names, d[()].ravel())}
    raise RuntimeError('event_table not found')

def read_log(path):
    wb = openpyxl.load_workbook(path, read_only=True, data_only=True)
    info = {}
    for sheet in ['agl_t_issue', 'agl_t']:          # agl_t overrides
        rows = list(wb[sheet].iter_rows(values_only=True)); h = rows[0]
        for r in rows[1:]:
            d = dict(zip(h, r))
            if not d.get('session') or d.get('probe_idx') is None: continue
            area = d.get('area_label_sec') or d.get('area_label')
            info[(d['session'], int(d['probe_idx']))] = dict(
                monkey=d['monkey'], area=str(area), spacing=d.get('electrode_spacing'),
                faulty=d.get('faulty_ch'), log_sheet=sheet)
    return info

def parse_faulty(x):
    if x is None: return []
    out = []
    for t in str(x).replace(';', ',').split(','):
        t = t.strip()
        try:
            v = int(float(t))
            if v > 1 or ',' in str(x): out.append(v)   # a lone '1' in Walt's log means "none"
        except ValueError: pass
    return out

# --------------------------------------------------------------------------- helpers
def patch(lfp, bad_local):
    lfp = lfp.copy(); n = lfp.shape[1]
    good = [c for c in range(n) if c not in bad_local]
    for c in bad_local:
        lo = max([g for g in good if g < c], default=None); hi = min([g for g in good if g > c], default=None)
        nb = [g for g in (lo, hi) if g is not None]
        lfp[:, c] = lfp[:, nb].mean(1)
    return lfp

def detect_bad(x):
    """Data-driven bad contacts: flat (SD < 20% of probe median), decoupled from both neighbours (broadband r < 0.5),
    or gamma-band (30-95 Hz) neighbour correlation < half the probe median or < median - 3 MAD."""
    seg = x[::5][:200000].astype(float)
    sd = seg.std(0); C = np.corrcoef(seg.T); n = x.shape[1]; bad = set(np.where(sd < 0.2 * np.median(sd))[0])
    for c in range(n):
        nb = [C[c, k] for k in (c - 1, c + 1) if 0 <= k < n and k not in bad]
        if nb and max(nb) < 0.5: bad.add(c)
    from scipy.signal import butter, sosfiltfilt
    yg = sosfiltfilt(butter(4, [30, 95], 'band', fs=FS, output='sos'), seg, axis=0)
    Cg = np.corrcoef(yg.T)
    nbg = np.array([np.mean([Cg[c, k] for k in (c - 1, c + 1) if 0 <= k < n]) for c in range(n)])
    med = np.median(nbg); mad = 1.4826 * np.median(np.abs(nbg - med)) + 1e-6
    bad |= set(np.where((nbg < 0.5 * med) | (nbg < med - 3 * mad))[0])
    return sorted(int(b) for b in bad)

def equalise_gain(x):
    """Remove contact gain/impedance differences: scale each contact to equal 1-8 Hz RMS (volume-conducted
    delta is near-uniform over the 2-3 mm probe span). Returns (x, gains)."""
    from scipy.signal import butter, sosfiltfilt
    sos = butter(3, [1, 8], 'band', fs=FS, output='sos')
    seg = x[:min(len(x), 300000)].astype(float)
    r = sosfiltfilt(sos, seg, axis=0).std(0)
    g = r / np.median(r)
    return x / g[None, :], g

def epochs(x, onsets, pre, post):
    on = onsets[np.isfinite(onsets)].astype(int)
    on = on[(on - pre >= 0) & (on + post < x.shape[0])]
    idx = on[:, None] + np.arange(-pre, post)[None, :]
    return x[idx]                                    # trials x time x ch

def band_mean(f, P, bands):
    m = np.zeros(len(f), bool)
    for lo, hi in bands: m |= (f >= lo) & (f <= hi)
    return P[:, m].mean(1)

def best_window_G(ab, g):
    n = len(ab); best = (0., None)
    for L in range(MIN_WIN, n + 1):
        for s in range(0, n - L + 1):
            x = np.arange(s, s + L)
            ra, rg = linregress(x, ab[x]), linregress(x, g[x])
            G = ra.rvalue**2 * rg.rvalue**2 if np.sign(ra.slope) != np.sign(rg.slope) else 0.
            if G > best[0]: best = (G, (s, s + L, ra, rg))
    return best

def csd(erp, h_um):
    # erp: time x ch ; Vaknin end-padding, 3-point Hamming spatial smoothing, units uV/mm^2 (relative)
    v = np.concatenate([erp[:, :1], erp, erp[:, -1:]], 1)
    w = np.array([0.23, 0.54, 0.23])
    vs = np.stack([np.convolve(v[t], w, 'same') for t in range(v.shape[0])])
    vs[:, 0], vs[:, -1] = v[:, 0], v[:, -1]
    h = h_um / 1000.
    return -(vs[:, 2:] - 2 * vs[:, 1:-1] + vs[:, :-2]) / h**2

# --------------------------------------------------------------------------- per-probe analysis
def block_index(x):
    """Gamma-band (30-95,105-145 Hz) inter-contact correlation; best two-block split.
    Returns (index, split_after_contact): within-block minus between-block mean correlation."""
    from scipy.signal import butter, sosfiltfilt
    seg = x[:min(len(x), 300000)].astype(float)                     # first 5 min
    sos1 = butter(4, [30, 95], 'band', fs=FS, output='sos'); sos2 = butter(4, [105, 145], 'band', fs=FS, output='sos')
    y = sosfiltfilt(sos1, seg, axis=0) + sosfiltfilt(sos2, seg, axis=0)
    C = np.corrcoef(y.T); n = len(C); best = (-1., None)
    for k in range(3, n - 2):
        A, B = slice(0, k), slice(k, n)
        w = np.r_[C[A, A][np.triu_indices(k, 1)], C[B, B][np.triu_indices(n - k, 1)]].mean()
        bi = w - C[A, B].mean()
        if bi > best[0]: best = (bi, k)
    return best[0], best[1], C

def analyse_probe(lfp16, onsets, spacing):
    res = {}
    # 1. spectrolaminar motif on the continuous recording (2-s Welch segments)
    f, P = welch(lfp16.astype(float), fs=FS, nperseg=2000, noverlap=1000, axis=0)   # f x ch
    P = P.T
    rel = P / P.max(0, keepdims=True)
    fm = (f >= 2) & (f <= 145) & ~((f > 47) & (f < 53)) & ~((f > 97) & (f < 103))
    L = np.log10(P[:, fm]); L -= L.mean(1, keepdims=True)
    shape_r = np.array([np.corrcoef(L[c], np.median(L[[k for k in range(max(0, c - 2), min(16, c + 3)) if k != c]], 0))[0, 1] for c in range(16)])
    res['n_odd_contacts'] = int(np.sum(shape_r < 0.9))
    ab, g = band_mean(f, rel, [AB_BAND]), band_mean(f, rel, G_BANDS)
    G, win = best_window_G(ab, g)
    null = np.array([best_window_G(ab[q], g[q])[0] for q in (rng.permutation(len(ab)) for _ in range(N_PERM))])
    res.update(G=G, p_motif=(1 + np.sum(null >= G)) / (1 + N_PERM), r_ab_gamma=float(np.corrcoef(ab, g)[0, 1]))
    if win:
        s, e, ra, rg = win
        xc = (rg.intercept - ra.intercept) / (ra.slope - rg.slope)
        res.update(win_start=s + 1, win_end=e, slope_ab=ra.slope, slope_g=rg.slope, crossover_ch=float(xc + 1),
                   crossover_inside=bool(s + 1 <= xc <= e - 2),
                   gamma_peak_end='low-index' if rg.slope < 0 else 'high-index')
    # 2. laminar block structure of gamma-band correlations
    bi, split, C = block_index(lfp16)
    res.update(block_index=float(bi), block_split_after=int(split))
    # 3. evoked CSD (secondary) on all element onsets, 0-150 ms
    on = onsets[np.isfinite(onsets)]
    el = (on[:, None] + np.array([0, 563, 1126, 1689, 2252])[None, :]).ravel()
    ep = epochs(lfp16, el, 100, 300).astype(float)
    ep -= ep[:, :100].mean(1, keepdims=True)
    erp = ep.mean(0); w = slice(100, 250)
    sem = ep.std(0).mean() / np.sqrt(len(ep))
    res['evoked_snr'] = float(np.abs(erp[w]).mean() / sem)
    res['min_erp_corr'] = float(np.nanmin(np.corrcoef(erp[w].T)))
    c_odd, c_even = csd(ep[0::2].mean(0), spacing), csd(ep[1::2].mean(0), spacing)
    res['csd_splithalf_r'] = float(np.corrcoef(c_odd[w].ravel(), c_even[w].ravel())[0, 1])
    res['n_events'] = int(len(ep))
    prof = dict(f=f, rel=rel, ab=ab, g=g, erp=erp, csd=csd(erp, spacing), gammaC=C, null=null)
    return res, prof

def classify(r, area):
    motif = r['G'] >= 0.3 and r['p_motif'] < 0.05 and r.get('crossover_inside', False)
    evoked_testable = r['evoked_snr'] >= 3
    evoked = evoked_testable and r['csd_splithalf_r'] >= 0.6
    block = (r['block_index'] >= 0.15 and 'crossover_ch' in r and
             abs((r['block_split_after'] + 0.5) - r['crossover_ch']) <= 2)
    r.update(motif_pass=bool(motif), evoked_testable=bool(evoked_testable), evoked_pass=bool(evoked), block_corroborates=bool(block))
    if motif and (evoked or block):   r['verdict'] = 'perpendicular'
    elif motif or evoked:             r['verdict'] = 'probable'
    else:                             r['verdict'] = 'no laminar signature'
    r['confidence'] = 'low' if r['n_odd_contacts'] >= 3 else 'normal'
    if str(area).lower().startswith(('hpc', 'hipp')): r['verdict'] = 'hippocampus (n/a)'
    return r

# --------------------------------------------------------------------------- main
def main():
    mat_dir, log_xlsx, out_dir = sys.argv[1:4]
    sessions = [x for x in sys.argv[4:] if x != '--force'] or sorted(x[:-4] for x in os.listdir(mat_dir) if x.endswith('.mat'))
    log = read_log(log_xlsx)
    os.makedirs(os.path.join(out_dir, 'probe_data'), exist_ok=True)
    for s in sessions:
        out_json = os.path.join(out_dir, 'probe_data', s + '.json')
        if os.path.exists(out_json) and '--force' not in sys.argv: continue
        try:
            with h5py.File(os.path.join(mat_dir, s + '.mat'), 'r') as f:
                ev = read_event_table(f)
                lfp = f['lfp'][()].astype(np.float32)
            onsets = ev['stimulusOnset_ms']
            rows = []
            for p in (1, 2):
                meta = log.get((s, p), dict(monkey=s.split('-')[0], area='unknown', spacing=150, faulty=None, log_sheet='none'))
                spacing = float(meta['spacing'] or 150)
                ch0 = (p - 1) * 16
                bad = detect_bad(lfp[:, ch0:ch0 + 16])
                x = patch(lfp[:, ch0:ch0 + 16], bad)
                x, gains = equalise_gain(x)
                r, prof = analyse_probe(x, onsets, spacing)
                r = classify(r, meta['area'])
                r.update(session=s, probe=p, monkey=meta['monkey'], area=meta['area'],
                         spacing_um=spacing, patched_ch=[b + 1 + ch0 for b in bad], max_gain_dev=float(np.max(np.abs(np.log2(gains)))), log_sheet=meta['log_sheet'])
                np.savez_compressed(os.path.join(out_dir, 'probe_data', f'{s}_p{p}.npz'), **prof)
                rows.append(r)
            json.dump(rows, open(out_json, 'w'), default=float)
            print(s, [(r['area'], round(r['G'], 2), round(r['p_motif'], 3), round(r['block_index'], 2), r['block_split_after'],
                       round(r['evoked_snr'], 1), r['verdict']) for r in rows], flush=True)
        except Exception as e:
            print(s, 'ERROR', repr(e), flush=True)

if __name__ == '__main__':
    main()
