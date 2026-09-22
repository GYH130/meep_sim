"""Record preservation and existing runtime, without solver execution or secrets."""
import hashlib
import json
from pathlib import Path
import platform
import subprocess

import ti2d  # Route runtime caches to the data disk before importing plot libraries.

import h5py
import matplotlib
import numpy
import scipy

from ti2d.common import atomic_json, digest_file, utcnow
from ti2d.queue import cgroup_resources

stage = Path(__file__).resolve().parent
if str(stage) != '/root/autodl-tmp/meep_sim/meep_screen':
    raise RuntimeError('Server only')
root = stage.parent
snapshots = sorted((stage / 'snapshots').glob('*/manifest.json'))
if not snapshots:
    raise RuntimeError('Snapshot required')
manifest_path = snapshots[0]
manifest = json.loads(manifest_path.read_text())
changed, missing = [], []
for relative, digest in manifest['files_sha256'].items():
    path = root / relative
    if not path.is_file():
        missing.append(relative)
    elif digest_file(path) != digest:
        changed.append(relative)
git_status = subprocess.check_output(['git', 'status', '--porcelain=v1', '-uall'], cwd=root, text=True)
legacy_dirty = [line for line in git_status.splitlines() if not line[3:].startswith('meep_screen/')]
state = {'created_utc': utcnow(), 'snapshot': str(manifest_path.parent),
         'legacy_files_checked': len(manifest['files_sha256']),
         'legacy_changed_since_snapshot': changed, 'legacy_missing_since_snapshot': missing,
         'legacy_git_status': legacy_dirty, 'not_found_at_start': manifest['not_found'],
         'runtime': {'python': platform.python_version(), 'numpy': numpy.__version__, 'scipy': scipy.__version__,
                     'h5py': h5py.__version__, 'matplotlib': matplotlib.__version__,
                     'environment': '/root/miniconda3/envs/meep_sim', 'installed_or_upgraded': False},
         'resources': cgroup_resources(), 'local_FDTD_executed': False}
atomic_json(stage / 'state_audit.json', state)
print(json.dumps({k: state[k] for k in ('legacy_files_checked', 'legacy_changed_since_snapshot',
                                      'legacy_missing_since_snapshot', 'not_found_at_start', 'resources')}))
if changed or missing:
    raise SystemExit(2)
