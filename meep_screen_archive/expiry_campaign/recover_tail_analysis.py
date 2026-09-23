"""Offline, evidence-preserving recovery of tail diagnostic reports.

No Meep import, solver launch, cache update, or qualification promotion. Original
measurements and reports are read-only; generated artifacts live exclusively in
reports_recovered/. An off-grid final observation is retained in the original
NPZ and explicitly excluded from FFT, never interpolated or silently clipped.
"""
from pathlib import Path
import argparse
import hashlib
import io
import json
import math
import os
import re

import numpy as np

from tail_analysis import replay_audit, spectrum
from ti2d.common import atomic_json, digest_file, utcnow


def regular_selection(times, fields, expected_dt, allow_terminal_partial=True):
    """Return audited row indices; reject any defect except one final short step."""
    t, values = np.asarray(times), np.asarray(fields)
    if (t.ndim != 1 or len(t) < 2 or np.iscomplexobj(t)
            or values.ndim < 1 or values.shape[0] != len(t)):
        raise ValueError('WAVEFORM_SHAPE_INVALID')
    if not np.isfinite(t).all() or not np.isfinite(values).all():
        raise ValueError('WAVEFORM_NONFINITE')
    dt = float(expected_dt)
    if not math.isfinite(dt) or dt <= 0:
        raise ValueError('EXPECTED_DT_INVALID')
    steps = np.diff(t)
    if np.any(steps <= 0):
        raise ValueError('TIME_DUPLICATE_OR_REVERSED')
    bad = np.flatnonzero(~np.isclose(steps, dt, atol=1e-8, rtol=1e-7))
    excluded = []
    keep = np.arange(len(t))
    if len(bad):
        if not (allow_terminal_partial and len(bad) == 1
                and bad[0] == len(steps)-1 and 0 < steps[-1] < dt):
            raise ValueError('INTERIOR_OR_UNSUPPORTED_TIME_GRID_DEFECT')
        keep = keep[:-1]
        excluded = [dict(index=len(t)-1, time=float(t[-1]),
                         preceding_time=float(t[-2]), actual_dt=float(steps[-1]),
                         reason='Final explicit observation lies after last regular sample but before next scheduled sample; excluded from FFT only; original retained')]
    return keep, dict(original_samples=len(t), regular_samples=len(keep),
                      expected_dt=dt, excluded_samples=excluded,
                      regular_interval=[float(t[keep[0]]), float(t[keep[-1]])],
                      resampled=False, interpolated=False, original_evidence_modified=False)


def checked_npz(payload, expected_sha256):
    if hashlib.sha256(payload).hexdigest() != expected_sha256:
        raise ValueError('WAVEFORM_CHECKSUM_FAILED')
    with np.load(io.BytesIO(payload), allow_pickle=False) as archive:
        arrays = {name: archive[name] for name in archive.files}
    for name, value in arrays.items():
        if not np.issubdtype(value.dtype, np.number) or not np.isfinite(value).all():
            raise ValueError('CHUNK_ARRAY_NONFINITE_OR_NONNUMERIC:'+name)
    return arrays


def resolve_status(worker=None,process=None):
    """Supervisor cause is authoritative; preserve the worker's own observation."""
    worker_status=worker.get('status') if worker is not None else None
    process_status=process.get('status') if process is not None else None
    if process is not None and not isinstance(process_status,str):
        raise ValueError('SUPERVISOR_STATUS_INVALID')
    if worker is not None and not isinstance(worker_status,str):
        raise ValueError('WORKER_STATUS_INVALID')
    return dict(status=process_status if process is not None else worker_status or 'NOT_RUN_OR_NO_RESULT',
                worker_status=worker_status,process_status=process_status,
                status_authority='process.json' if process is not None else 'result.json' if worker is not None else 'no process or worker result')


def load_waveforms(read_bytes, arm, manifest, metadata, state):
    ts, values, line_ts, line_values = [], [], [], []
    xline, seen = None, set()
    for item in manifest:
        name = item['file']
        if not re.fullmatch(r'waveforms_\d{4,}\.npz', name) or name in seen:
            raise ValueError('MANIFEST_FILENAME_INVALID_OR_DUPLICATE')
        seen.add(name)
        z = checked_npz(read_bytes(arm+'/'+name), item['sha256'])
        t, v = z['point_t'], z['point_fields']
        if (t.ndim != 1 or not len(t)
                or v.shape != (len(t), len(metadata['probe_names']), len(metadata['point_component_names']))):
            raise ValueError('CHUNK_POINT_SHAPE_MISMATCH')
        if float(t[0]) != item['first_time'] or float(t[-1]) != item['last_time']:
            raise ValueError('MANIFEST_TIME_MISMATCH')
        ts.append(t); values.append(v)
        lt = z['line_t']
        if lt.ndim != 1:
            raise ValueError('CHUNK_LINE_TIME_SHAPE_MISMATCH')
        if len(lt):
            x, lv = z['x_back'], z['back_Ey']
            if x.ndim != 1 or lv.shape != (len(lt),len(x)):
                raise ValueError('CHUNK_LINE_SHAPE_MISMATCH')
            if xline is not None and not np.array_equal(xline,x):
                raise ValueError('LINE_COORDINATES_CHANGED')
            if not np.isin(lt,t).all():
                raise ValueError('LINE_SAMPLE_NOT_IN_POINT_GRID')
            xline = x
            line_ts.append(lt); line_values.append(lv)
    if not ts:
        raise ValueError('NO_WAVEFORM_CHUNKS')
    t, v = np.concatenate(ts), np.concatenate(values)
    keep, audit = regular_selection(t,v,metadata['point_dt'])
    if not np.isclose(t[-1],state['simulation_time'],atol=1e-8,rtol=0):
        # New probe versions may put the terminal observation in a separate
        # checksummed artifact, leaving all waveform chunks strictly uniform.
        gap=float(state['simulation_time'])-float(t[-1])
        if not (0 < gap < float(metadata['point_dt'])):
            raise ValueError('FINAL_SAMPLE_RESULT_TIME_MISMATCH')
        terminal_meta=json.loads(read_bytes(arm+'/terminal_sample.json'))
        if terminal_meta['file'] != 'terminal_sample.npz':
            raise ValueError('TERMINAL_SAMPLE_FILENAME_INVALID')
        terminal=checked_npz(read_bytes(arm+'/terminal_sample.npz'),terminal_meta['sha256'])
        if (terminal['point_t'].shape != (1,)
                or terminal['point_fields'].shape != (1,)+v.shape[1:]
                or float(terminal['point_t'][0]) != float(state['simulation_time'])
                or float(terminal_meta['time']) != float(state['simulation_time'])):
            raise ValueError('TERMINAL_SAMPLE_RESULT_MISMATCH')
        audit['separate_terminal_sample']=dict(time=float(terminal['point_t'][0]),
            gap_after_regular_sample=gap,sha256=terminal_meta['sha256'],
            excluded_from_fft=True,reason='Explicit off-grid final observation stored separately; original retained')
    audit['all_chunk_checksums_verified'] = True
    audit['chunks_verified'] = len(manifest)
    audit['all_numeric_arrays_finite'] = True
    lines = None
    if line_ts:
        lt, lv = np.concatenate(line_ts), np.concatenate(line_values)
        # A final point could also trigger a line capture. Apply the same explicit
        # terminal-only rule; every retained line must lie on the regular grid.
        lk, la = regular_selection(lt,lv,metadata['line_dt'])
        period=metadata['case']['geometry']['period_um']
        if (len(xline)<2 or not np.all(np.diff(xline)>0)
                or not np.allclose(np.diff(xline),period/len(xline),atol=1e-8,rtol=1e-7)):
            raise ValueError('SPATIAL_GRID_NOT_ONE_REGULAR_PERIOD')
        audit['line_sampling'] = la
        lines = (lt[lk],lv[lk],xline)
    return t[keep],v[keep],lines,audit


def boundary_comparison(arms, waveforms, control_reproduced, floor):
    """Compare completed arms on identical retained timestamps, never interpolate."""
    names=('pml6_control','pml12_test')
    if not control_reproduced or any(arms.get(n,{}).get('status')!='DIAGNOSTIC_COMPLETE'
                                    or arms.get(n,{}).get('worker_status')!='DIAGNOSTIC_COMPLETE' for n in names):
        return None
    if any(n not in waveforms for n in names):
        raise ValueError('COMPARISON_WAVEFORMS_MISSING')
    if not math.isfinite(floor) or floor<=0:
        raise ValueError('COMPARISON_DENOMINATOR_INVALID')
    end=min(float(waveforms[n]['point_t'][-1]) for n in names)
    begin=max(600.,end-200.,max(float(waveforms[n]['source_end']) for n in names))
    selected={}
    for name in names:
        w=waveforms[name]
        t=np.asarray(w['point_t']);values=np.asarray(w['point_fields'])
        keep,_=regular_selection(t,values,w['point_dt'],allow_terminal_partial=False)
        use=keep[(t[keep]>=begin)&(t[keep]<=end)]
        if len(use)<128:
            raise ValueError('COMPARISON_POST_SOURCE_WINDOW_INSUFFICIENT')
        selected[name]=(t[use],values[use])
    tc,vc=selected[names[0]];tt,vt=selected[names[1]]
    if not np.array_equal(tc,tt):
        raise ValueError('COMPARISON_TIME_GRIDS_MISMATCH')
    stats={};spectra={}
    for name in names:
        t,values=selected[name];w=waveforms[name]
        index=w['probe_names'].index('back_air')
        ratios=np.sum(np.abs(values[:,index,:2])**2,axis=1)/floor
        stats[name]=dict(mean_intensity_ratio=float(ratios.mean()),
                         max_intensity_ratio=float(ratios.max()))
        spectra[name]={}
        for probe in ('front_air','back_air'):
            i=w['probe_names'].index(probe)
            _,_,detail=spectrum(t,values[:,i,:2])
            spectra[name][probe]=detail
    ctrl,test=stats[names[0]],stats[names[1]]
    return dict(status='MATCHED_DIAGNOSTIC_COMPARISON',qualification=False,
                requested_interval=[begin,end],interval=[float(tc[0]),float(tc[-1])],
                samples=len(tc),sample_dt=float(tc[1]-tc[0]),
                matching='Exact identical retained timestamps; no interpolation or resampling',
                normalization='(|Ex|²+|Ey|²) / H reference stop denominator',
                control_mean_intensity_ratio=ctrl['mean_intensity_ratio'],
                thick_pml_mean_intensity_ratio=test['mean_intensity_ratio'],
                thick_over_control=test['mean_intensity_ratio']/ctrl['mean_intensity_ratio'] if ctrl['mean_intensity_ratio']>0 else None,
                control_max_intensity_ratio=ctrl['max_intensity_ratio'],
                thick_pml_max_intensity_ratio=test['max_intensity_ratio'],
                thick_over_control_max=test['max_intensity_ratio']/ctrl['max_intensity_ratio'] if ctrl['max_intensity_ratio']>0 else None,
                matched_window_spectra=spectra,
                interpretation='PML-thickness sensitivity diagnostic only. Differences do not by themselves establish causality, convergence, or optical qualification.')


def analyze_evidence(read_bytes):
    """Analyze from a read-only byte provider, usable with directory or tar data."""
    evidence_hashes = {}

    def read_json(name):
        raw=read_bytes(name)
        evidence_hashes[name]=hashlib.sha256(raw).hexdigest()
        return json.loads(raw)

    h=read_json('inputs/H_result.json')
    htrace=read_json('inputs/H_stop.json')['trace']
    original_status=read_json('status.json')
    floor=float(h['stop']['normalization']['denominator_floor'])
    if not math.isfinite(floor) or floor<=0:
        raise ValueError('REFERENCE_DENOMINATOR_INVALID')
    result=dict(schema_version=1, generated_utc=utcnow(), qualification=False,
                recovery_status='OFFLINE_DIAGNOSTIC_RECOVERED',
                original_run_status=original_status, H_status=h['status'],
                control_reproduced=False, arms={}, evidence_sha256=evidence_hashes,
                scope='Specified Ordal Route A Ti model; 2D slanted slot; oxide=false; diagnostic only',
                limitations=['No new FDTD or qualification.',
                             'A recovered spectrum does not satisfy convergence or replay coverage.',
                             'A PML-thickness comparison is a sensitivity diagnostic, not a causal or qualification proof.',
                             'Signed complex field spectra and spatial amplitudes are not diffraction powers.'])
    plotting={'source_end':float(h['source_end_time']), 'H_trace':htrace,
              'denominator':floor,'arms':{}}
    waveforms={}
    for arm in ('pml6_control','pml12_test'):
        try:
            process=read_json(arm+'/process.json')
        except FileNotFoundError:
            process=None
        try:
            state=read_json(arm+'/result.json')
        except FileNotFoundError:
            result['arms'][arm]=dict(resolve_status(None,process),qualification=False)
            continue
        metadata=read_json(arm+'/metadata.json')
        trace=read_json(arm+'/power_trace.json')
        manifest=read_json(arm+'/waveform_manifest.json')
        t,fields,lines,audit=load_waveforms(read_bytes,arm,manifest,metadata,state)
        if metadata['point_component_names'][:2] != ['Ex','Ey']:
            raise ValueError('E_COMPONENT_ORDER_MISMATCH')
        a=dict(resolve_status(state,process),qualification=False,simulation_time=state['simulation_time'],
               elapsed_s=state['elapsed_s'],sampling_audit=audit,spectra={})
        if arm=='pml6_control':
            result['control_replay']=replay_audit(htrace,trace,floor)
            result['control_reproduced']=result['control_replay']['passed']
            result['control_replay']['required_coverage_time']=1000.
        tail=t>=max(600.,float(metadata['source_end'])+75.)
        plot_arm={'trace':trace,'spectra':{}}
        if tail.sum()>=128:
            for probe in ('front_air','back_air'):
                i=metadata['probe_names'].index(probe)
                freqs,power,detail=spectrum(t[tail],fields[tail,i,:2])
                a['spectra'][probe]=detail
                plot_arm['spectra'][probe]=(freqs,power)
        else:
            a['spectral_status']='INSUFFICIENT_POST_SOURCE_SAMPLES'
        if lines is not None:
            lt,lv,x=lines
            use=lt>=max(600.,float(metadata['source_end'])+75.)
            if use.sum()>1:
                case=metadata['case']; period=case['geometry']['period_um']
                kx=-math.sin(math.radians(case['emission_theta_deg']))/case['wavelength_um']
                modes=np.fft.fft(lv[use]*np.exp(-2j*np.pi*kx*x)[None,:],axis=1)/len(x)
                weights=np.mean(np.abs(modes)**2,axis=0)
                ms=np.fft.fftfreq(len(x))*len(x)
                ix=np.argsort(weights)[::-1][:6]
                a['back_Ey_spatial_orders']=[dict(m=int(round(ms[i])),
                    mean_squared_amplitude=float(weights[i]),
                    rayleigh_cutoff=float(abs(kx+ms[i]/period))) for i in ix]
                a['spatial_interval']=[float(lt[use][0]),float(lt[use][-1])]
                a['spatial_order_note']='Ey complex amplitudes after Bloch phase removal; not powers or shares'
        result['arms'][arm]=a
        plotting['arms'][arm]=plot_arm
        waveforms[arm]=dict(point_t=t,point_fields=fields,point_dt=metadata['point_dt'],
                           probe_names=metadata['probe_names'],source_end=metadata['source_end'])
    comparison=boundary_comparison(result['arms'],waveforms,result['control_reproduced'],floor)
    if comparison is not None:
        result['boundary_comparison']=comparison
    else:
        result['boundary_comparison_status']='NOT_EVALUATED: requires completed control and PML12 arms plus original control reproduction gate'
    return result,plotting


def write_outputs(out,result,plotting):
    os.environ['MPLCONFIGDIR']=str(Path(__file__).resolve().parent/'cache/matplotlib')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    out.mkdir(exist_ok=True)
    fig,axes=plt.subplots(1,2,figsize=(12,4.4),constrained_layout=True)
    source_end=plotting['source_end']; floor=plotting['denominator']
    hrows=[r for r in plotting['H_trace'] if r['time']>=source_end]
    for name in ('front_air','back_air','cavity_0.5'):
        axes[0].semilogy([r['time'] for r in hrows],[r['probe_ratio'][name] for r in hrows],label=name)
    axes[0].set_title('H: original unqualified post-source tail')
    for arm,p in plotting['arms'].items():
        rows=[r for r in p['trace'] if r['time']>=source_end]
        if rows:
            axes[1].semilogy([r['time'] for r in rows],[r['probe_power']['back_air']/floor for r in rows],label=arm)
    axes[1].set_title('Measured back-side tail (diagnostic)')
    for ax in axes:
        ax.set(xlabel='Simulation time (um/c)',ylabel='Intensity / reference stop denominator')
        ax.axhline(1e-4,color='black',ls='--',label='stop threshold')
        ax.grid(alpha=.2); ax.legend(fontsize=8)
    fig.savefig(out/'tail_comparison.png',dpi=160);plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(12,4.4),constrained_layout=True)
    for probe,ax in zip(('front_air','back_air'),axes):
        for arm,p in plotting['arms'].items():
            if probe not in p['spectra']:continue
            f,power=p['spectra'][probe]; ix=np.argsort(f)
            ax.semilogy(f[ix],power[ix]/max(float(power.max()),1e-300),label=arm)
        ax.set(title=probe+' regular-segment spectrum',xlabel='Signed frequency (1/um)',
               ylabel='Relative Ex/Ey spectral intensity',xlim=(-.5,.5),ylim=(1e-10,2))
        ax.grid(alpha=.2)
        if ax.lines:ax.legend(fontsize=8)
    fig.savefig(out/'tail_spectra.png',dpi=160);plt.close(fig)
    lines=['# Tail diagnostic offline report recovery','',
           'Diagnostic only. Original measurements, original status and original reports are unchanged.',
           'Recovered FFT uses verified regular samples; no interpolation or resampling.','',
           'Original run state: '+result['original_run_status']['state']+'.',
           'Original H qualification: '+result['H_status']+'.',
           'Control reproduction gate passed: '+str(result['control_reproduced'])+'.','']
    for name,a in result['arms'].items():
        lines.append('- '+name+': '+a['status']+'; qualification=false.')
        if a.get('process_status') is not None:
            lines.append('  Authoritative supervisor status=%s; worker status=%s.'%(a['process_status'],a.get('worker_status')))
        if 'sampling_audit' in a:
            s=a['sampling_audit']
            lines.append('  Regular samples: %s; original samples: %s.'%(s['regular_samples'],s['original_samples']))
            for excluded in s['excluded_samples']:
                lines.append('  FFT-only exclusion: terminal t=%s, interval=%s instead of %s; original NPZ retained.'%(excluded['time'],excluded['actual_dt'],s['expected_dt']))
            if 'separate_terminal_sample' in s:
                terminal=s['separate_terminal_sample']
                lines.append('  Separate terminal observation t=%s retained in original checksummed NPZ; excluded from FFT.'%terminal['time'])
    if 'control_replay' in result:
        r=result['control_replay']
        lines+=['','Matched H through t=%s; required t>=1000; maximum normalized error=%s.'%(r['through_time'],r['max_power_error_over_reference'])]
    if 'boundary_comparison' in result:
        b=result['boundary_comparison']
        lines+=['','Exact matched PML6/PML12 interval: %s; %s samples.'%(b['interval'],b['samples']),
                'Mean back intensity / H reference: control=%s; thick PML=%s; thick/control=%s.'%(b['control_mean_intensity_ratio'],b['thick_pml_mean_intensity_ratio'],b['thick_over_control']),
                'Maximum back intensity / H reference: control=%s; thick PML=%s; thick/control=%s.'%(b['control_max_intensity_ratio'],b['thick_pml_max_intensity_ratio'],b['thick_over_control_max']),
                'Matched-window signed complex spectra are recorded in diagnostic_summary.json.',
                b['interpretation']]
    else:
        lines+=['',result['boundary_comparison_status']]
    lines+=['','These diagnostics do not establish a PML-boundary cause or a qualified optical size trend.',
            'Details and evidence hashes: diagnostic_summary.json.','']
    (out/'handoff.md').write_text('\n'.join(lines))
    atomic_json(out/'diagnostic_summary.json',result)


def recover(run,output=None):
    run=Path(run).resolve()
    if not run.is_dir():raise ValueError('RUN_DIRECTORY_REQUIRED')

    def read_bytes(name):
        path=(run/name).resolve()
        if run not in path.parents:raise ValueError('EVIDENCE_PATH_OUTSIDE_RUN')
        return path.read_bytes()

    out=Path(output).resolve() if output is not None else run/'reports_recovered'
    if run not in out.parents or not out.relative_to(run).parts[0].startswith('reports_recovered'):
        raise ValueError('OUTPUT_MUST_BE_SEPARATE_DIRECTORY_INSIDE_RUN')
    try:
        result,plotting=analyze_evidence(read_bytes)
        result['recovery_implementation_sha256']=digest_file(__file__)
        write_outputs(out,result,plotting)
    except Exception as exc:
        out.mkdir(exist_ok=True)
        atomic_json(out/'recovery_error.json',dict(status='RECOVERY_REJECTED',
                    error=str(exc),qualification=False,original_evidence_modified=False))
        raise
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('run',type=Path)
    parser.add_argument('--output',type=Path,default=None)
    args=parser.parse_args()
    recovered=recover(args.run,args.output)
    print(json.dumps(dict(status=recovered['recovery_status'],qualification=False,
                          control_reproduced=recovered['control_reproduced'],
                          reports=str(args.output.resolve() if args.output else args.run.resolve()/'reports_recovered'))))
