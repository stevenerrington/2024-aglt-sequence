"""Apply final perpendicularity rules to per-probe metrics from laminar_perpendicularity.py.
Usage: python3 summarise_perpendicularity.py <results_dir>
Writes <results_dir>/perpendicularity_summary.csv (one row per probe) and updates each probe's verdict in the JSONs.

Rules (per 16-contact probe)
  motif     : G >= 0.3, no artefact boundary at the crossover, permutation p < 0.05, alpha-beta/gamma crossover inside the fitted window, AND gamma
              peaking at the same probe end as the monkey's majority (contact order differs between rigs:
              Troy/Plexon gamma at low index, Walt/TDT at high index). Opposite orientation -> 'inverted motif'.
  evoked    : element-onset evoked SNR >= 3, CSD split-half r >= 0.6, AND evoked-LFP polarity reversal across
              depth (min inter-contact ERP correlation < -0.3).
  block     : gamma-band correlation block split within 2 contacts of the spectral crossover (block index >= 0.15).
  verdict   : perpendicular = motif & (evoked | block); probable = motif or evoked; else no laminar signature.
  artefacts : a gamma-correlation discontinuity (min adjacent r < 0.6 x probe median) within 2 contacts of the
              crossover marks an artefact/out-of-brain boundary -> motif rejected; block structure not counted on
              probes with any discontinuity; evoked CSD with > 40% of its energy on one contact rejected.
  NOTE      : evoked metrics use uncorrected stimulusOnset_ms (no audio latency); Troy's evoked metrics are
              uninformative (no audio channel, trial-by-trial latency unknown).
  confidence: low if >= 2 contacts have a spectral shape deviating from neighbours (RMS log-dev > 0.1) or the
              probe has a gamma-correlation discontinuity.
"""
import sys, json, glob, numpy as np, pandas as pd
res = sys.argv[1]
AREA = {'AntAud': 'auditory', 'PostAud(R)': 'auditory', 'R': 'auditory', 'R/A1': 'auditory', '44.0': 'frontal',
        '45.0': 'frontal', '46.0': 'frontal', 'FOP': 'frontal', 'Hipp CA1?': 'hpc', 'dSTS?': 'other'}
files = sorted(glob.glob(f'{res}/probe_data/*.json'))
R = [r for f in files for r in json.load(open(f))]
for r in R:
    z = np.load(f"{res}/probe_data/{r['session']}_p{r['probe']}.npz"); f, rel = z['f'], z['rel']
    fm = (f >= 2) & (f <= 145) & ~((f > 47) & (f < 53)) & ~((f > 97) & (f < 103))
    L = np.log10(rel[:, fm] + 1e-12); L -= L.mean(1, keepdims=True)
    dev = [np.sqrt(np.mean((L[c] - np.median(L[[k for k in range(max(0, c - 2), min(16, c + 3)) if k != c]], 0)) ** 2)) for c in range(16)]
    r['odd_contacts'] = [i + 1 for i, v in enumerate(dev) if v > 0.1]
    adj = np.diag(z['gammaC'], 1)
    r['disc_ratio'] = float(adj.min() / np.median(adj)); r['disc_after'] = int(adj.argmin()) + 1
    e = (z['csd'][100:250] ** 2).sum(0)
    r['csd_concentration'] = float(np.nanmax(e) / np.nansum(e)) if np.nansum(e) > 0 else np.nan
    r['area_group'] = AREA.get(r['area'], r['area'])
d = pd.DataFrame(R)
d['boundary_at_crossover'] = (d.disc_ratio < 0.6) & ((d.disc_after + 0.5 - d.crossover_ch).abs() <= 2)
raw_motif = (d.G >= 0.3) & (d.p_motif < 0.05) & d.crossover_inside.astype(bool) & ~d.boundary_at_crossover
ref = d[raw_motif].groupby('monkey').gamma_peak_end.agg(lambda s: s.value_counts().idxmax()).to_dict()
print('reference gamma end per monkey:', ref, d[raw_motif].groupby('monkey').gamma_peak_end.value_counts().to_dict())
for r in R:
    r['boundary_at_crossover'] = bool(r['disc_ratio'] < 0.6 and abs(r['disc_after'] + .5 - r.get('crossover_ch', np.nan)) <= 2)
    raw = r['G'] >= 0.3 and r['p_motif'] < 0.05 and r.get('crossover_inside', False) and not r['boundary_at_crossover']
    r['motif_orientation_ok'] = bool(raw and r.get('gamma_peak_end') == ref.get(r['monkey']))
    motif = r['motif_orientation_ok']
    evoked = r['evoked_snr'] >= 3 and r['csd_splithalf_r'] >= 0.6 and r['min_erp_corr'] < -0.3 and r['csd_concentration'] <= 0.4
    block = r['disc_ratio'] >= 0.6 and r['block_index'] >= 0.15 and np.isfinite(r.get('crossover_ch', np.nan)) and abs(r['block_split_after'] + .5 - r['crossover_ch']) <= 2
    r.update(motif_pass=bool(motif), evoked_pass=bool(evoked), block_corroborates=bool(block))
    if r['area_group'] == 'hpc':                 v = 'hippocampus (n/a)'
    elif motif and (evoked or block):            v = 'perpendicular'
    elif motif or evoked:                        v = 'probable'
    elif raw:                                    v = 'inverted motif (check)'
    else:                                        v = 'no laminar signature'
    r['verdict'] = v
    r['confidence'] = 'low' if (len(r['odd_contacts']) >= 2 or r['disc_ratio'] < 0.6) else 'normal'
for f in files:
    s = json.load(open(f))[0]['session']
    json.dump([r for r in R if r['session'] == s], open(f, 'w'), default=lambda o: o if not isinstance(o, np.generic) else o.item())
d = pd.DataFrame(R)
cols = ['session', 'probe', 'monkey', 'area', 'area_group', 'verdict', 'confidence', 'motif_pass', 'evoked_pass',
        'block_corroborates', 'G', 'p_motif', 'crossover_ch', 'gamma_peak_end', 'block_index', 'block_split_after',
        'evoked_snr', 'csd_splithalf_r', 'min_erp_corr', 'csd_concentration', 'disc_ratio', 'disc_after', 'boundary_at_crossover', 'odd_contacts', 'patched_ch', 'n_events', 'log_sheet']
d[cols].round(3).to_csv(f'{res}/perpendicularity_summary.csv', index=False)
print(pd.crosstab([d.monkey, d.area_group], d.verdict, margins=True))
print(pd.crosstab(d.verdict, d.confidence))
