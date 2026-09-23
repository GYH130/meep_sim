"""Offline waveform audit and plots. Never imports Meep or changes qualification."""
from pathlib import Path
import argparse
import json
import math
import numpy as np
from ti2d.common import atomic_json, digest_file


def replay_audit(h_trace, replay, floor):
    old = {round(x['time'], 5): x for x in h_trace}
    errors, tail_errors, times = [], [], []
    for row in replay:
        ref = old.get(round(row['time'], 5))
        if ref is None:
            continue
        diffs = [abs(row['probe_power'][n]-v)/floor for n,v in ref['probe_power'].items()]
        errors += diffs; times.append(row['time'])
        if row['time'] >= 600:
            tail_errors += diffs
    return dict(compared_samples=len(times), through_time=max(times, default=0),
        max_power_error_over_reference=max(errors, default=None),
        max_tail_error_over_reference=max(tail_errors, default=None),
        passed=bool(tail_errors and max(times)>=1000 and max(errors)<1e-4 and max(tail_errors)<1e-5),
        scope='Same source, material, grid and PML; passive-monitor omission replay only')


def spectrum(t, fields):
    t, fields = np.asarray(t), np.asarray(fields)
    if len(t)<128 or fields.shape[0]!=len(t):
        raise ValueError('INSUFFICIENT_WAVEFORM')
    dt = float(np.median(np.diff(t)))
    if not np.allclose(np.diff(t), dt, atol=1e-8, rtol=1e-7):
        raise ValueError('NONUNIFORM_TIME_GRID')
    # Keep DC; a quasistatic residual is one possible explanation to examine.
    window = np.hanning(len(t))
    a = np.fft.fft(fields * window[:, None], axis=0) / window.sum()
    power = np.sum(np.abs(a)**2, axis=1)
    freq = np.fft.fftfreq(len(t), dt)
    selected = []
    for i in np.argsort(power)[::-1]:
        if all(abs(freq[i]-freq[j]) > 3/(len(t)*dt) for j in selected):
            selected.append(int(i))
        if len(selected)==6:
            break
    return freq, power, dict(sample_dt=dt, frequency_resolution=1/(len(t)*dt),
        nyquist=1/(2*dt), interval=[float(t[0]),float(t[-1])],
        peaks=[dict(signed_frequency=float(freq[i]), absolute_frequency=float(abs(freq[i])),
                    relative_power=float(power[i]/max(float(power.max()),1e-300))) for i in selected],
        convention='Complex fields; retain both signed frequencies; Hann window; Ex/Ey power sum; no DC subtraction')


def analyze(root, run):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    root, run = Path(root), Path(run)
    h = json.loads((run/'inputs/H_result.json').read_text())
    htrace = json.loads((run/'inputs/H_stop.json').read_text())['trace']
    floor = h['stop']['normalization']['denominator_floor']
    plan = json.loads((run/'plan.json').read_text())
    result = dict(qualification=False, H_status=h['status'], H_last_time=htrace[-1]['time'],
                  H_last_probe_ratio=htrace[-1]['probe_ratio'], arms={}, control_reproduced=False)
    figs = run/'reports'; figs.mkdir(exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.3), constrained_layout=True)
    for name in ('front_air','back_air','cavity_0.5'):
        axes[0].semilogy([r['time'] for r in htrace], [r['probe_ratio'][name] for r in htrace], label=name)
    axes[0].axhline(1e-4, color='black', ls='--', label='qualification threshold')
    axes[0].axvline(h['source_end_time'], color='gray', ls=':')
    axes[0].set(xlabel='Simulation time (um/c)', ylabel='Intensity / stop denominator', title='Previous H: unqualified tail')
    axes[0].legend(fontsize=8)
    spectra_plot, spectrum_axes = plt.subplots(1, 2, figsize=(12,4.3), constrained_layout=True)
    curves = {}
    for arm in ('pml6_control','pml12_test'):
        path = run/arm
        if not (path/'result.json').exists():
            result['arms'][arm] = dict(status='NOT_COMPLETE'); continue
        state = json.loads((path/'result.json').read_text())
        metadata = json.loads((path/'metadata.json').read_text())
        trace = json.loads((path/'power_trace.json').read_text())
        arm_result = dict(status=state['status'], simulation_time=state['simulation_time'],
                          elapsed_s=state['elapsed_s'], diagnostic_only=True)
        if arm=='pml6_control':
            result['control_replay'] = replay_audit(htrace, trace, floor)
            result['control_reproduced'] = result['control_replay']['passed']
        tt = np.array([r['time'] for r in trace])
        yy = np.array([r['probe_power']['back_air']/floor for r in trace])
        curves[arm] = tt, yy
        axes[1].semilogy(tt, yy, label=arm)
        manifests = json.loads((path/'waveform_manifest.json').read_text())
        time_chunks, value_chunks = [], []
        line_t, line_chunks = [], []
        xline = None
        for item in manifests:
            p = path/item['file']
            if digest_file(p)!=item['sha256']:
                raise ValueError('WAVEFORM_CHECKSUM_FAILED')
            with np.load(p,allow_pickle=False) as z:
                time_chunks.append(z['point_t']); value_chunks.append(z['point_fields'])
                if len(z['line_t']):
                    line_t.append(z['line_t'])
                    line_chunks.append(z['back_Ey'])
                    xline=z['x_back']
        ts, val = np.concatenate(time_chunks), np.concatenate(value_chunks)
        tail = ts>=max(600,metadata['source_end']+75)
        if tail.sum()>=128:
            arm_result['spectra']={}
            for pn, ax in zip(('front_air','back_air'),spectrum_axes):
                i = metadata['probe_names'].index(pn)
                freqs, powr, detail = spectrum(ts[tail],val[tail,i,:2])
                arm_result['spectra'][pn]=detail
                order=np.argsort(freqs)
                ax.semilogy(freqs[order],powr[order]/max(powr.max(),1e-300),label=arm)
                ax.set(xlim=(-.5,.5),ylim=(1e-10,2),xlabel='Signed frequency (1/um)',ylabel='Relative spectral power',title=pn+' tail spectrum')
                ax.axvline(-1/10.5,color='gray',ls=':'); ax.legend(fontsize=8)
        if line_t:
            lt, lv = np.concatenate(line_t), np.concatenate(line_chunks)
            use=lt>=600
            kx=-math.sin(math.radians(30))/10.5
            if use.sum()>1:
                modes=np.fft.fft(lv[use]*np.exp(-2j*np.pi*kx*xline)[None,:],axis=1)/len(xline)
                weights=np.mean(np.abs(modes)**2,axis=0)
                m=np.fft.fftfreq(len(xline))*len(xline)
                selected=np.argsort(weights)[::-1][:6]
                arm_result['back_Ey_spatial_orders']=[dict(m=int(round(m[i])),
                    mean_squared_amplitude=float(weights[i]), rayleigh_cutoff=float(abs(kx+m[i]/50))) for i in selected]
                arm_result['spatial_order_note']='Field amplitudes, not diffraction powers; periodic phase removed; no qualification inference'
        result['arms'][arm]=arm_result
    axes[1].axhline(1e-4,color='black',ls='--')
    axes[1].set(xlabel='Simulation time (um/c)',ylabel='Back |E|² / H reference peak',title='Boundary-thickness diagnostic')
    if curves: axes[1].legend(fontsize=8)
    fig.savefig(figs/'tail_comparison.png',dpi=150); plt.close(fig)
    spectra_plot.savefig(figs/'tail_spectra.png',dpi=150); plt.close(spectra_plot)
    if all(a in curves for a in ('pml6_control','pml12_test')) and result['control_reproduced']:
        a,b=curves['pml6_control'],curves['pml12_test']
        end=min(a[0][-1],b[0][-1]); begin=max(600.,end-200.)
        use=(a[0]>=begin)&(a[0]<=end)
        if use.sum()>10:
            ctrl=float(np.mean(a[1][use])); test=float(np.mean(np.interp(a[0][use],b[0],b[1])))
            result['boundary_comparison']=dict(interval=[begin,float(end)],control_mean_intensity_ratio=ctrl,
                thick_pml_mean_intensity_ratio=test,thick_over_control=test/ctrl if ctrl>0 else None,
                note='Sensitivity to PML thickness supports boundary involvement; it does not by itself qualify an optical result')
    atomic_json(figs/'diagnostic_summary.json',result)
    lines=['# I: slanted-slot tail diagnostic','',
           'This is a diagnostic, not a qualified size study. Frozen Ordal Route A; 2D; oxide=false.',
           '',f"H remains {h['status']}; last simulation time {htrace[-1]['time']:.3f}.",
           f"PML6 passive replay reproduced: {result['control_reproduced']}.",'']
    for arm, data in result['arms'].items():
        lines.append(f"- {arm}: {data['status']}; t={data.get('simulation_time')}; wall_s={data.get('elapsed_s')}")
    if 'boundary_comparison' in result:
        lines += ['', 'Boundary comparison: '+json.dumps(result['boundary_comparison'])]
    lines += ['', 'Full spectra and spatial-order evidence: diagnostic_summary.json.',
              'Figures: tail_comparison.png, tail_spectra.png.',
              'Required next decision: distinguish boundary-sensitive tail from a persistent material/geometry mode using waveform evidence; no automatic follow-up solver.',
              'No source bandwidth, material, geometry, Courant or stop threshold has been relaxed. No full timestep-state restart exists.',
              'AI/API calls in background queue: zero. Actual AI tokens and rental RMB cost unavailable; no invented estimate.','']
    (figs/'handoff.md').write_text('\n'.join(lines))
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(); p.add_argument('run'); a=p.parse_args()
    analyze(Path(__file__).resolve().parent,Path(a.run))
