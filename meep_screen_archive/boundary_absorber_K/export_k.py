"""Checkpoint K evidence without network access, Meep execution or legacy edits.

Call only between owned solver tasks, after their processes have exited. Each
generation is immutable; export_ready.json is published last and atomically.
"""
from datetime import datetime, timezone
import argparse
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import tarfile
import tempfile

ROOT = Path(__file__).resolve().parent
PREFIX = Path('meep_screen_archive/boundary_absorber_K')
TEXT_SUFFIXES = {'.py', '.md', '.json', '.jsonl', '.csv', '.txt', '.yml', '.yaml', '.log', '.patch', '.diff', '.toml'}
COMPACT_SUFFIXES = {'.py', '.md', '.json', '.csv', '.txt', '.yml', '.yaml', '.png', '.svg', '.toml'}
SECRETS = [rb'github_pat_[A-Za-z0-9_]{25,}', rb'gh[pousr]_[A-Za-z0-9]{25,}',
           rb'sk-(?:ant-|proj-)[A-Za-z0-9_-]{20,}',
           rb'-----BEGIN (?:OPENSSH |RSA |EC )?PRIVATE KEY-----',
           rb'https?://[^\s/:]+:[^\s/@]+@']


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def atomic_json(path, obj):
    path = Path(path)
    temporary = path.with_name('.' + path.name + '.' + str(os.getpid()))
    with temporary.open('w') as stream:
        json.dump(obj, stream, indent=2, ensure_ascii=False, allow_nan=False)
        stream.write('\n')
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(temporary, path)


def scan(path):
    # Selected binary NPZ/HDF5 evidence is never put in the public archive.
    if path.suffix not in TEXT_SUFFIXES | {'.svg'}:
        return
    overlap = b''
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            data = overlap + block
            if any(re.search(pattern, data) for pattern in SECRETS):
                raise ValueError('SECRET_SCAN_BLOCKED:' + path.name)
            overlap = data[-1024:]


def selected_files(root):
    paths = []
    for folder in ('ti2d', 'tests', 'assets', 'inputs', 'runs', 'cache/references'):
        base = root / folder
        if base.exists():
            paths.extend(base.rglob('*'))
    paths.extend(root.glob('*'))
    selected = {}
    for path in sorted(set(paths)):
        relative = path.relative_to(root)
        if any(part.startswith('.') or part == '__pycache__' for part in relative.parts):
            continue
        if path.is_symlink():
            raise ValueError('SYMLINK_NOT_ARCHIVED:' + str(relative))
        if not path.is_file() or path.suffix in {'.pyc', '.lock'}:
            continue
        if len(relative.parts) == 1 and path.suffix not in COMPACT_SUFFIXES | {'.log'}:
            continue
        scan(path)
        selected[str(relative)] = path
    return selected


def verify_results(files):
    verified = []
    for relative, path in files.items():
        if path.name != 'result.json':
            continue
        obj = json.loads(path.read_text())
        status = obj.get('status', 'UNKNOWN')
        if status == 'QUALIFIED':
            evidence = obj.get('evidence_sha256')
            if not evidence:
                raise ValueError('QUALIFIED_WITHOUT_EVIDENCE:' + relative)
            for name, expected in evidence.items():
                evidence_path = path.parent / name
                if (evidence_path.is_symlink()
                        or not evidence_path.resolve().is_relative_to(path.parent.resolve())
                        or not evidence_path.is_file() or sha(evidence_path) != expected):
                    raise ValueError('QUALIFIED_EVIDENCE_HASH_FAILED:' + relative)
        verified.append({'path': relative, 'status': status})
    return verified


def main(final=False, label=None):
    root = ROOT
    exports = root / 'exports'
    exports.mkdir(parents=True, exist_ok=True)
    with (exports / 'export.lock').open('a+') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        now = datetime.now(timezone.utc)
        generation = now.strftime('%Y%m%dT%H%M%S%fZ')
        label = label or ('final' if final else 'checkpoint')
        if not re.fullmatch(r'[A-Za-z0-9_-]{1,64}', label):
            raise ValueError('UNSAFE_CHECKPOINT_LABEL')
        files = selected_files(root)
        verification = verify_results(files)
        fingerprints = {name: sha(path) for name, path in files.items()}
        raw_name = 'K_raw_evidence_' + generation + '.tar.gz'
        compact_name = 'K_compact_' + generation + '.tar.gz'
        with tempfile.TemporaryDirectory(prefix='snapshot_', dir=exports) as temporary:
            temporary = Path(temporary)
            metadata = {
                'generation': generation, 'generated_utc': now.isoformat(),
                'label': label, 'final': bool(final), 'case_results': verification,
                'scope': 'Individual case evidence only; no automatic production qualification.',
                'solver_stop_cst': '2026-09-24T20:30:00+08:00',
                'provider_expiry_cst_user_confirmed_approximate': '2026-09-24T21:53:00+08:00',
                'raw_files_sha256': fingerprints,
            }
            ledger = root.parent / 'budget_ledger.json'
            if ledger.is_file():
                # Freeze the shared ledger bytes; never change the active ledger.
                (temporary / 'budget_ledger_snapshot.json').write_bytes(ledger.read_bytes())
            atomic_json(temporary / 'K_SNAPSHOT.json', metadata)
            raw_temp = temporary / raw_name
            with tarfile.open(raw_temp, 'w:gz') as archive:
                for name, path in files.items():
                    archive.add(path, arcname='boundary_absorber_K/' + name, recursive=False)
                archive.add(temporary / 'K_SNAPSHOT.json', arcname='boundary_absorber_K/K_SNAPSHOT.json')
                if (temporary / 'budget_ledger_snapshot.json').exists():
                    archive.add(temporary / 'budget_ledger_snapshot.json', arcname='boundary_absorber_K/budget_ledger_snapshot.json')
            # A changing source cannot be published as a coherent checkpoint.
            if any(sha(path) != fingerprints[name] for name, path in files.items()):
                raise RuntimeError('EVIDENCE_CHANGED_DURING_EXPORT')
            dest = temporary / PREFIX
            dest.mkdir(parents=True)
            inventory = []

            def add(path, name):
                if path.stat().st_size > 8 * 1024**2:
                    raise ValueError('COMPACT_FILE_TOO_LARGE:' + name)
                scan(path)
                target = dest / name
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(path, target)
                inventory.append({'path': name, 'sha256': sha(target), 'bytes': target.stat().st_size})

            for name, path in files.items():
                if ('cache' not in Path(name).parts and path.suffix in COMPACT_SUFFIXES
                        and not path.name.startswith('heartbeat')):
                    add(path, name)
            add(temporary / 'K_SNAPSHOT.json', 'K_SNAPSHOT.json')
            if (temporary / 'budget_ledger_snapshot.json').exists():
                add(temporary / 'budget_ledger_snapshot.json', 'budget_ledger_snapshot.json')
            atomic_json(dest / 'K_FILES.json', inventory)
            compact_temp = temporary / compact_name
            with tarfile.open(compact_temp, 'w:gz') as archive:
                archive.add(dest, arcname=str(PREFIX))
            raw_sha, compact_sha = sha(raw_temp), sha(compact_temp)
            os.replace(raw_temp, exports / raw_name)
            os.replace(compact_temp, exports / compact_name)
            shutil.copyfile(exports / compact_name, temporary / 'compact_latest.tar.gz')
            os.replace(temporary / 'compact_latest.tar.gz', exports / 'compact_latest.tar.gz')
            ready = {'generation': generation, 'final': bool(final), 'label': label,
                     'raw_file': raw_name, 'raw_sha256': raw_sha,
                     'compact_file': compact_name, 'compact_sha256': compact_sha,
                     'generated_utc': now.isoformat(), 'case_results': verification}
            atomic_json(exports / ('export_' + generation + '.json'), ready)
            atomic_json(exports / 'export_ready.json', ready)
            return ready


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--final', action='store_true')
    parser.add_argument('--label')
    args = parser.parse_args()
    print(json.dumps(main(final=args.final, label=args.label), indent=2))
