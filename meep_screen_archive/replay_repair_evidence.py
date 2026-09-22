"""Replay saved stop traces without importing Meep or changing their status."""
import json
from pathlib import Path
import sys

from ti2d.common import atomic_json, digest_file
from ti2d.diagnostics import StopTracker

stage = Path(__file__).resolve().parent
if sys.platform != 'linux' or str(stage) != '/root/autodl-tmp/meep_sim/meep_screen':
    raise RuntimeError('SERVER_EVIDENCE_REQUIRED')
old = stage / 'runs/stage_20260921_B/tasks/flat_r16_p_plus30/attempt_1'
reference = json.loads((old / 'reference_stop.json').read_text())
structure = json.loads((old / 'structure_stop.json').read_text())
source_end = structure['source_end']
source_peak = max(v['envelope_power'] for v in reference['trace'] if v['time'] <= source_end)
records = []
for threshold in (1e-4, 1e-6):
    tracker = StopTracker(source_end, window=20, consecutive_windows=3,
                          decay_threshold=threshold, rt_tolerance=1e-3,
                          reference_intensity_peak=source_peak,
                          reference_provenance='saved runB reference source-on non-PML |E|² peak',
                          reference_floor_fraction=1.0)
    for sample in structure['trace']:
        tracker.sample(sample['time'], sample['probe_power'], sample['envelope_power'], sample['R'], sample['T'])
        if tracker.converged:
            break
    records.append({'threshold': threshold, 'replayed_stop': tracker.summary(),
                    'physical_status': 'UNQUALIFIED_ORIGINAL_REFLECTION_IS_INVALID'})
result = {'status': 'OFFLINE_REPLAY_ONLY', 'old_status': structure['stop_reason'],
          'old_last_time': structure['last_time'], 'source_end': source_end,
          'reference_excitation_peak': source_peak, 'replays': records,
          'evidence_sha256': {p.name: digest_file(p) for p in (old / 'reference_stop.json', old / 'structure_stop.json')},
          'warning': 'Neither replay nor changed stop logic qualifies the old run. No original file modified.'}
path = stage / 'repair_diagnostics/stop_replay_repair1.json'
if path.exists():
    raise FileExistsError('Keep the previous replay instead of overwriting')
atomic_json(path, result)
print(json.dumps({'path': str(path), 'old_last_time': result['old_last_time'],
                  'new_first_stop_times': [r['replayed_stop']['last_time'] for r in records],
                  'converged': [r['replayed_stop']['converged'] for r in records],
                  'physical_qualification': False}))
