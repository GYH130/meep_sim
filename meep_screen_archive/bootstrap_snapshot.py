"""Read-only legacy snapshot; all writes remain under the new stage directory."""
import hashlib
import json
import os
import subprocess
import tarfile
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path('/root/autodl-tmp/meep_sim')
STAGE = ROOT / 'meep_screen'


def run(*args):
    return subprocess.check_output(args, cwd=ROOT)


def main():
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    dest = STAGE / 'snapshots' / stamp
    dest.mkdir(parents=True, exist_ok=False)
    (dest / 'git_status.txt').write_bytes(run('git', 'status', '--porcelain=v1', '-uall'))
    (dest / 'git_head.txt').write_bytes(run('git', 'rev-parse', 'HEAD'))
    (dest / 'git_branch.txt').write_bytes(run('git', 'branch', '--show-current'))
    (dest / 'unstaged.patch').write_bytes(run('git', 'diff', '--binary'))
    (dest / 'staged.patch').write_bytes(run('git', 'diff', '--cached', '--binary'))
    candidates = ['src', 'scripts', 'tests', 'configs', 'reports', 'results', 'review', 'logs', 'sources']
    found, missing, hashes = [], [], {}
    for name in candidates:
        p = ROOT / name
        (found if p.exists() else missing).append(name)
    extra = [p for p in ROOT.iterdir() if p.is_file() and (p.suffix in ('.md', '.txt', '.yml', '.yaml', '.toml') or p.name == '.gitignore')]
    selected = [ROOT / n for n in found] + extra
    for target in selected:
        for p in ([target] if target.is_file() else sorted(target.rglob('*'))):
            if p.is_file() and not p.is_symlink() and '__pycache__' not in p.parts and '.ipynb_checkpoints' not in p.parts:
                hashes[str(p.relative_to(ROOT))] = hashlib.sha256(p.read_bytes()).hexdigest()
    with tarfile.open(dest / 'legacy_evidence.tar.gz', 'w:gz') as archive:
        for name in hashes:
            archive.add(ROOT / name, arcname=name, recursive=False)
    info = {'created_utc': stamp, 'found': found, 'not_found': missing, 'files_sha256': hashes,
            'scope': 'Existing project code, configs, reports, results and logs; no .git, credentials or home configuration.'}
    (dest / 'manifest.json').write_text(json.dumps(info, indent=2, ensure_ascii=False) + '\n')
    print(json.dumps({'snapshot': str(dest), 'file_count': len(hashes), 'not_found': missing,
                      'archive_bytes': (dest / 'legacy_evidence.tar.gz').stat().st_size}))


if __name__ == '__main__':
    main()
