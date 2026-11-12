"""Manifest-driven server M0 geometry runs with host/WSL resource monitoring.

No SSH, scheduler mutation, package installation, or source-model edits.
An existing output directory is never resumed or overwritten implicitly.
"""
from __future__ import annotations
import argparse
from collections import defaultdict, deque
import concurrent.futures as cf
import datetime as dt
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import re
import signal
import shutil
import subprocess
import threading
import time
from resource_guard import Guard, stop_owned

ROOT=Path(__file__).resolve().parent
KIND='IDEAL_GEOMETRY_PILOT_NOT_TA2_SPECTRAL_OR_FULL_WAVE_VALIDATION'
DIMENSIONS = ('height_um', 'bottom_width_um', 'upper_width_um', 'pitch_um')

def utc():return dt.datetime.now(dt.timezone.utc).isoformat()
def digest(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,v, *, compact=False):
    p=Path(p);temp=p.with_suffix(p.suffix+'.tmp')
    temp.write_text(json.dumps(v,ensure_ascii=False,indent=None if compact else 2,
                               separators=(',', ':') if compact else None,allow_nan=False)+'\n')
    os.replace(temp,p)

def integer(value, label, minimum, maximum):
    if type(value) is not int or not minimum <= value <= maximum:
        raise ValueError(f'{label} must be an integer in [{minimum}, {maximum}]')
    return value


def number(value, label):
    if type(value) not in (int, float) or not math.isfinite(value):
        raise ValueError(f'{label} must be finite and numeric')
    return float(value) or 0.0


def derived_count(record, key, value, label):
    if key in record and (type(record[key]) is not int or record[key] != value):
        raise ValueError(f'{label}.{key} is {record[key]!r}, expected {value}')
    record[key] = value


def angular_axis(phase, key, limit, label):
    values = phase.get(key)
    if not isinstance(values, list) or not values:
        raise ValueError(f'{label}.{key} must be a nonempty angle list')
    angles = [number(value, f'{label}.{key}') for value in values]
    if len(set(angles)) != len(angles) or any(not 0 <= angle < limit for angle in angles):
        raise ValueError(f'{label}.{key} must contain distinct angles in [0, {limit})')
    return sorted(angles)


def manifest():
    m = json.loads((ROOT / 'manifest.json').read_text())
    if m['schema'] != 2 or m['kind'] != KIND:
        raise ValueError('Unsupported M0 batch manifest')
    if m['engine_source'] != 'ray_material.cpp' or digest(ROOT / 'ray_material.cpp') != m['engine_source_sha256']:
        raise ValueError('Ray engine source digest mismatch')
    integer(m['max_bounces'], 'max_bounces', 1, 256)
    integer(m['seed_base'], 'seed_base', 0, (1 << 63) - 1)
    if not isinstance(m['shapes'], list) or not m['shapes']:
        raise ValueError('shapes must be a nonempty list')
    shapes_by_id = {}
    used_ids = set()
    for shape in m['shapes']:
        shape_id = shape['id']
        if not isinstance(shape_id, str) or not re.fullmatch(r'[A-Za-z][A-Za-z0-9_-]{0,47}', shape_id):
            raise ValueError(f'Unsafe geometry ID: {shape_id!r}')
        if shape_id.casefold() in used_ids:
            raise ValueError(f'Duplicate geometry ID: {shape_id}')
        if shape['family'] not in ('A0', 'A3'):
            raise ValueError(f'Unsupported family for {shape_id}')
        h, b, d, p = (number(shape[key], f'{shape_id}.{key}') for key in DIMENSIONS)
        if not (h > 0 and 0 < b < d < p):
            raise ValueError(f'{shape_id} requires H > 0 and 0 < b < D < P')
        integer(shape['curve_segments'], f'{shape_id}.curve_segments', 1, 4096)
        shape.update(zip(DIMENSIONS, (h, b, d, p)))
        shape.update(b_over_D=b / d, H_over_D=h / d, D_over_P=d / p)
        shapes_by_id[shape_id] = shape
        used_ids.add(shape_id.casefold())
    derived_count(m, 'shape_count', len(shapes_by_id), 'manifest')

    if not isinstance(m['phases'], dict) or not m['phases']:
        raise ValueError('phases must be a nonempty object')
    total = 0
    for name, phase in m['phases'].items():
        if not re.fullmatch(r'[a-z][a-z0-9_]{0,31}', name):
            raise ValueError(f'Unsafe phase name: {name!r}')
        label = f'phases.{name}'
        ids = phase.get('shape_ids')
        if not isinstance(ids, list) or not ids or any(not isinstance(s, str) or s not in shapes_by_id for s in ids):
            raise ValueError(f'{label}.shape_ids must list known geometry IDs')
        if len(set(ids)) != len(ids):
            raise ValueError(f'Duplicate geometry in {label}.shape_ids')
        theta = angular_axis(phase, 'theta_deg', 90, label)
        phi = angular_axis(phase, 'phi_deg', 360, label)
        phase.update(theta_deg=theta, phi_deg=phi)
        # Polar azimuth is undefined at normal incidence, so evaluate it once.
        count = sum(1 if t == 0 else len(phi) for t in theta)
        integer(phase['rays_per_condition'], f'{label}.rays_per_condition', 1, 1_000_000_000)
        derived_count(phase, 'shape_count', len(ids), label)
        derived_count(phase, 'angles_per_shape', count, label)
        derived_count(phase, 'condition_count', len(ids) * count, label)
        total += phase['condition_count']
    derived_count(m, 'condition_count', total, 'manifest')
    return m


def angle_token(value, width):
    # repr of a float round-trips exactly. No periods in paths used with with_suffix.
    text = str(int(value)) if value.is_integer() else repr(value)
    integer_part, point, remainder = text.partition('.')
    return (integer_part.zfill(width) + (f'p{remainder}' if point else '')).replace('+', '')


def jobs(m, phase_name):
    phase = m['phases'][phase_name]
    shapes = {s['id']: s for s in m['shapes']}
    for shape_id in phase['shape_ids']:
        shape = shapes[shape_id]
        dims = [shape[key] for key in DIMENSIONS]
        for theta in phase['theta_deg']:
            for phi in ([0.0] if theta == 0 else phase['phi_deg']):
                material = json.dumps([m['seed_base'], dims, theta, phi], separators=(',', ':')).encode()
                seed = int.from_bytes(hashlib.sha256(material).digest()[:8], 'big') & ((1 << 63) - 1)
                yield dict(
                    id=f'{shape_id}_th{angle_token(theta, 2)}_ph{angle_token(phi, 3)}',
                    shape=shape,
                    angle=dict(theta_deg=theta, phi_deg=phi),
                    rays=phase['rays_per_condition'],
                    max_bounces=m['max_bounces'],
                    seed=seed,
                )

def command(j):
    s,a=j['shape'],j['angle']
    return [str(ROOT/'ray_material'),'--structure',s['family'],'--surface','M0',
            '--height-um',str(s['height_um']),'--bottom-width-um',str(s['bottom_width_um']),
            '--upper-width-um',str(s['upper_width_um']),'--pitch-um',str(s['pitch_um']),
            '--theta',str(a['theta_deg']),'--phi',str(a['phi_deg']),
            '--rays',str(j['rays']),'--segments',str(s['curve_segments']),
            '--max-bounces',str(j['max_bounces']),'--seed',str(j['seed'])]

def verify_result(j,v):
    s,a=j['shape'],j['angle']
    assert v['surface']=='M0' and v['structure']==s['family']
    for key in ('height_um','bottom_width_um','upper_width_um','pitch_um'):
        assert abs(v[key]-s[key])<1e-9
    assert v['theta_deg']==a['theta_deg'] and v['phi_deg']==a['phi_deg']
    assert v['rays']==j['rays'] and v['seed']==j['seed']
    assert v['curve_segments']==s['curve_segments'] and v['max_bounces']==j['max_bounces']
    assert v['mesh_validation'].startswith('PASS_')
    assert 0<=v['max_ray_closure_error']<1e-9
    for pol in ('s','p','unpolarized'):
        c=v['channels'][pol]
        for k in ('A_total','A_active','A_wall','escape','residual'):
            assert math.isfinite(c[k]) and -1e-10<=c[k]<=1+1e-10
        assert abs(c['A_total']+c['escape']+c['residual']-1)<1e-9
        assert abs(c['A_wall'])<1e-10 and abs(c['A_total']-c['A_active'])<1e-10
    assert abs(v['channels']['s']['A_total']-v['channels']['p']['A_total'])<1e-10
    for k in ('cap_terminated_rays','geometry_unresolved_rays','cutoff_terminated_rays'):
        assert isinstance(v[k],int) and 0<=v[k]<=j['rays']
    c=v['channels']['unpolarized']
    return dict(scope=KIND,energy_accounting='PASS',geometry_binding='PASS',
                quality='EXPLORATORY_UNRESOLVED_WARNING' if c['residual']>0 or v['geometry_unresolved_rays']>0 else 'IDEAL_SAMPLE_ACCOUNTED',
                absorption_interval=[c['A_total'],c['A_total']+c['residual']],
                residual=c['residual'],se_spatial_only=c['se']['A_total'],
                cap_terminated_rays=v['cap_terminated_rays'],geometry_unresolved_rays=v['geometry_unresolved_rays'])

def verify_build(m):
    if platform.system()!='Linux' or platform.machine() not in ('x86_64','AMD64'):
        raise RuntimeError('Run on server Linux/WSL x86_64, not Mac')
    b=json.loads((ROOT/'build_receipt.json').read_text())
    assert b['self_test']['self_test']=='PASS'
    assert b['source_sha256']==m['engine_source_sha256']==digest(ROOT/'ray_material.cpp')
    assert b['manifest_sha256']==digest(ROOT/'manifest.json')
    assert b['binary_sha256']==digest(ROOT/'ray_material')
    return b

def check_platform():
    if platform.system()!='Linux' or platform.machine() not in ('x86_64','AMD64'):
        raise RuntimeError('Server Linux/WSL x86_64 required')
    # Preserve the user requirement that scratch and output stay on D:.
    configured_root = os.environ.get('MEEP_RAY_STAGE_ROOT')
    if not configured_root or not ROOT.is_relative_to(Path(configured_root).resolve()):
        raise RuntimeError('Set MEEP_RAY_STAGE_ROOT to the designated compute staging parent')

def build(args,m):
    check_platform()
    if (ROOT/'ray_material').exists() or (ROOT/'build_receipt.json').exists():
        raise RuntimeError('Existing binary or receipt: inspect without overwriting')
    compiler=shutil.which(args.compiler)
    if not compiler:raise RuntimeError('Existing compiler unavailable; do not install system packages implicitly')
    scratch=ROOT/'scratch';scratch.mkdir(exist_ok=True)
    env=dict(os.environ,TMPDIR=str(scratch),TMP=str(scratch),TEMP=str(scratch),OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
    guard=Guard();records=[]
    def step(cmd,label):
        while True:
            admission=guard.sample([],os.getpid());records.append(dict(step=label,admission=admission));write(ROOT/'build_resource_evidence.json',records)
            write(ROOT/'build_status.json',dict(state='STARTING_'+label.upper() if admission['admit_new'] else 'WAITING_RESOURCES',
                phase=label,utc=utc(),controller_pid=os.getpid(),resources=admission))
            if admission['admit_new']:break
            if (ROOT/'STOP_NEW_DISPATCH').exists():raise RuntimeError('Stop requested while waiting for build resources')
            time.sleep(30)
        with (ROOT/(label+'.stdout')).open('x') as out,(ROOT/(label+'.stderr')).open('x') as err:
            p=subprocess.Popen(cmd,stdout=out,stderr=err,env=env,cwd=ROOT,start_new_session=True)
            last=time.monotonic()
            while p.poll() is None:
                if time.monotonic()-last>=5:
                    s=guard.sample([p],os.getpid());records.append(dict(step=label,runtime=s));write(ROOT/'build_resource_evidence.json',records);last=time.monotonic()
                    if s['stop_owned_reasons']:
                        stop_owned([p]);raise RuntimeError('Resource guard stopped owned '+label)
                time.sleep(.1)
            if p.returncode:raise RuntimeError(label+' failed; retained stdout/stderr')
    step([compiler,'-O3','-std=c++17',str(ROOT/'ray_material.cpp'),'-o',str(ROOT/'ray_material')],'compile')
    step([str(ROOT/'ray_material'),'--self-test'],'self_test')
    t=json.loads((ROOT/'self_test.stdout').read_text());assert t['self_test']=='PASS'
    with (ROOT/'ray_material').open('rb') as f:header=f.read(20)
    assert header[:4]==b'\x7fELF' and header[4:6]==b'\x02\x01' and int.from_bytes(header[18:20],'little')==62
    v=dict(kind='LINUX_X86_64_PILOT_BUILD',utc=utc(),platform=platform.platform(),machine=platform.machine(),
           source_sha256=digest(ROOT/'ray_material.cpp'),binary_sha256=digest(ROOT/'ray_material'),manifest_sha256=digest(ROOT/'manifest.json'),
           self_test=t,compiler=compiler,flags=['-O3','-std=c++17'],resource_evidence='build_resource_evidence.json')
    assert v['source_sha256']==m['engine_source_sha256']
    write(ROOT/'build_receipt.json',v)
    write(ROOT/'build_status.json',dict(state='COMPLETE_BUILD_SELFTEST',utc=utc(),source_sha256=v['source_sha256']))
    print(json.dumps({'build':'PASS','machine':v['machine'],'source_sha256':v['source_sha256']}))
    return 0

def run(args,m):
    check_platform()
    build=verify_build(m)
    phase = m['phases'][args.phase]
    out = Path(args.output).resolve()
    if ROOT not in out.parents:
        raise ValueError('Output must be inside the server D: stage')
    out.mkdir(parents=True, exist_ok=False)
    source = jobs(m, args.phase)
    pending = deque([next(source)] if args.action == 'smoke' else source)
    shape_ids = {job['shape']['id'] for job in pending}
    for shape_id in shape_ids:
        (out / 'conditions' / shape_id).mkdir(parents=True)
    args.workers=min(args.workers,len(os.sched_getaffinity(0)))
    scratch=ROOT/'scratch';scratch.mkdir(exist_ok=True)
    guard=Guard();snapshot=guard.sample([],os.getpid())
    total=len(pending);began=time.monotonic();started=utc();active={};owned={};lock=threading.Lock();stop=[]
    for signum in (signal.SIGTERM,signal.SIGINT):signal.signal(signum,lambda s,f:stop.append(s))
    completed_count = 0
    unresolved_warnings = 0
    failed = []
    by_shape_angle = defaultdict(lambda: defaultdict(list))
    worker_seconds = defaultdict(float)
    write(out/'run_manifest.json',dict(kind=KIND,manifest=m,build=build,manifest_sha256=digest(ROOT/'manifest.json'),
          driver_sha256=digest(Path(__file__)),run_action=args.action,phase=args.phase,condition_count=total,
          workers=args.workers,guard='Embedded host/WSL sampling; 100GiB start, 64GiB runtime, 112/128GiB owned RSS, two consecutive violations, CPU75%, healthy D >=100GiB, WSL available >=8GiB',
          constraints={'OMP_NUM_THREADS':1,'OPENBLAS_NUM_THREADS':1,'MKL_NUM_THREADS':1}))
    def one(j):
        base=out/'conditions'/j['shape']['id']/j['id'];stamp=utc();t0=time.monotonic();cmd=command(j)
        write(base.with_suffix('.request.json'),dict(job=j,command=cmd,source_sha256=m['engine_source_sha256'],started_utc=stamp))
        env=dict(os.environ,TMPDIR=str(scratch),TMP=str(scratch),TEMP=str(scratch),OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
        with base.with_suffix('.stdout').open('x') as stdout,base.with_suffix('.stderr').open('x') as stderr:
            p=subprocess.Popen(cmd,stdout=stdout,stderr=stderr,env=env,start_new_session=True)
            with lock:
                active[j['id']]=dict(pid=p.pid,started_utc=stamp,command=cmd)
                owned[j['id']]=p
            code=p.wait()
        try:
            if code:raise RuntimeError(f'engine exit {code}; see retained stderr')
            v=json.loads(base.with_suffix('.stdout').read_text());check=verify_result(j,v)
            record=dict(kind=KIND,job=j,source_sha256=m['engine_source_sha256'],binary_sha256=build['binary_sha256'],
                        started_utc=stamp,finished_utc=utc(),elapsed_seconds=time.monotonic()-t0,
                        exit_code=code,result=v,checks=check)
            write(base.with_suffix('.json'),record)
            return (j['angle']['phi_deg'], v['channels']['unpolarized']['A_total'],
                    check['residual'], check['se_spatial_only']), record['elapsed_seconds'], check['quality']
        finally:
            with lock:
                active.pop(j['id'],None);owned.pop(j['id'],None)
    def state(state_name):
        with lock:current=dict(active)
        write(out/'status.json',dict(kind=KIND,phase=args.phase,state=state_name,started_utc=started,heartbeat_utc=utc(),
              controller_pid=os.getpid(),condition_count=total,completed_conditions=completed_count,
              planned_geometry_count=phase['shape_count'],angles_per_shape=phase['angles_per_shape'],
              rays_per_condition=phase['rays_per_condition'],
              failed=failed,not_dispatched=len(pending),active=current,workers=args.workers,
              elapsed_seconds=time.monotonic()-began,stop_requested=bool(stop),resource_snapshot=snapshot,formal_qualification=False))
    state('STARTING')
    with cf.ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures={};last=0.0;last_resource=time.monotonic()
        with (out/'resource_samples.jsonl').open('a') as f:f.write(json.dumps(snapshot)+'\n')
        while pending or futures:
            if not stop and (out/'STOP_NEW_DISPATCH').exists():stop.append('STOP_NEW_DISPATCH')
            if time.monotonic()-last_resource>=(5 if futures else 30):
                with lock:handles=list(owned.values())
                snapshot=guard.sample(handles,os.getpid());last_resource=time.monotonic()
                with (out/'resource_samples.jsonl').open('a') as f:f.write(json.dumps(snapshot)+'\n')
                if snapshot['stop_owned_reasons']:
                    stop.append({'resource_reasons':snapshot['stop_owned_reasons'],'stopped_owned_pids':stop_owned(handles)})
            while pending and len(futures)<args.workers and not stop and not failed and snapshot['admit_new']:
                j=pending.popleft();futures[pool.submit(one,j)]=j
            if not futures:
                if stop or failed:break
                state('WAITING_RESOURCES');time.sleep(1);continue
            done,_=cf.wait(futures,timeout=1,return_when=cf.FIRST_COMPLETED)
            for f in done:
                j=futures.pop(f)
                try:
                    sample, elapsed, quality = f.result()
                    shape_id = j['shape']['id']
                    by_shape_angle[shape_id][j['angle']['theta_deg']].append(sample)
                    worker_seconds[shape_id] += elapsed
                    completed_count += 1
                    if quality == 'EXPLORATORY_UNRESOLVED_WARNING':
                        unresolved_warnings += 1
                except Exception as e:failed.append(dict(id=j['id'],error=repr(e),utc=utc()))
            if time.monotonic()-last>=2 or not futures:
                state('DRAINING' if failed or stop else 'RUNNING');last=time.monotonic()
    summaries = []
    expected = phase['angles_per_shape']
    for shape in m['shapes']:
        shape_id = shape['id']
        if shape_id not in by_shape_angle:
            continue
        angles = []
        for theta, samples in sorted(by_shape_angle[shape_id].items()):
            samples.sort()
            values = [sample[1] for sample in samples]
            residuals = [sample[2] for sample in samples]
            angles.append(dict(theta_deg=theta, azimuth_count=len(samples), mean_A=math.fsum(values)/len(samples),
                  min_A=min(values), max_A=max(values), max_residual=max(residuals),
                  mean_A_upper_bound=math.fsum(v+r for v,r in zip(values,residuals))/len(samples),
                  azimuth_samples=samples))
        summaries.append(dict(shape=shape, completed_angles=sum(a['azimuth_count'] for a in angles),
               expected_angles=expected, summed_worker_seconds=worker_seconds[shape_id], angles=angles))
    final_state='FAILED_PARTIAL' if failed else ('STOPPED_PARTIAL' if completed_count<total else 'COMPLETE_IDEAL_PILOT')
    summary=dict(kind=KIND,phase=args.phase,state=final_state,started_utc=started,finished_utc=utc(),
          elapsed_seconds=time.monotonic()-began,
          manifest_sha256=build['manifest_sha256'],engine_source_sha256=build['source_sha256'],
          binary_sha256=build['binary_sha256'],
          condition_count=total,completed_conditions=completed_count,
          completed_geometries=sum(s['completed_angles']==expected for s in summaries),
          planned_geometry_count=phase['shape_count'],angles_per_shape=expected,
          rays_per_condition=phase['rays_per_condition'],failed=failed,not_dispatched=len(pending),
          stop_reasons=stop, azimuth_sample_columns=['phi_deg','A_total','residual','se_spatial_only'],
          max_residual=max((a['max_residual'] for s in summaries for a in s['angles']),default=None),
          unresolved_warning_conditions=unresolved_warnings,
          geometry_summaries=summaries,formal_qualification=False,
          limitations=['Finite angular grid and spatial sampling; refine geometry ranking where differences are small',
                       'Residual is an explicit absorption interval, never counted as escape',
                       'Ideal black floor and specular mirror walls, with no spectral or temperature model'])
    write(out/'summary.json',summary,compact=True);state(final_state)
    print(json.dumps({k:summary[k] for k in ('phase','state','elapsed_seconds','completed_conditions','completed_geometries','max_residual')}))
    return 0 if final_state == 'COMPLETE_IDEAL_PILOT' else 1

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('action',choices=['validate-plan','commands','build','smoke','run','resources'])
    p.add_argument('--output',default=str(ROOT/'results'))
    p.add_argument('--workers',type=int,default=16)
    p.add_argument('--compiler',default='g++')
    p.add_argument('--phase',help='Phase name in manifest.json; defaults to the first listed phase')
    args=p.parse_args()
    if not 1<=args.workers<=16:p.error('workers must be 1..16')
    m=manifest()
    args.phase = args.phase or next(iter(m['phases']))
    if args.phase not in m['phases']:
        p.error(f'Unknown phase {args.phase!r}; choose from {", ".join(m["phases"])}')
    if args.action=='validate-plan':
        phase_counts = {name: {key: phase[key] for key in ('shape_count','angles_per_shape','condition_count','rays_per_condition')}
                        for name, phase in m['phases'].items()}
        print(json.dumps({'plan':'PASS','shapes':m['shape_count'],'conditions':m['condition_count'],
                          'phases':phase_counts,'source_sha256':m['engine_source_sha256']}))
    elif args.action=='commands':print(json.dumps([dict(job=j,command=command(j)) for j in jobs(m,args.phase)],indent=2))
    elif args.action=='resources':
        check_platform();print(json.dumps(Guard().sample([],os.getpid()),indent=2))
    elif args.action=='build':raise SystemExit(build(args,m))
    else:raise SystemExit(run(args,m))
