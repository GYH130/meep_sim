"""Generate an explicit, secret-checked upload inventory; never git-add or push."""
import argparse
import json
import re
from pathlib import Path

from ti2d.common import atomic_json, digest_file, utcnow


def main(run):
    stage = Path(__file__).resolve().parent
    run = Path(run).resolve()
    if stage not in run.parents:
        raise ValueError('Run must be inside this stage')
    allowed = [stage / 'README.md', stage / 'export_handoff.py', stage / 'solver_smoke.py']
    for directory, pattern in ((stage / 'ti2d', '*.py'), (stage / 'tests', '*.py'),
                               (stage / 'assets/material', '*'), (run / 'cases', '*.json'),
                               (run / 'reports', '*.md'), (run / 'reports', '*.json'),
                               (run / 'reports', '*.csv'), (run / 'reports/figures', '*.png')):
        allowed.extend(sorted(directory.glob(pattern)))
    allowed += [run / name for name in ('manifest.json', 'geometries.json', 'material_recheck.json')]
    # Entire credentials, not variable names or placeholder prefixes.
    patterns = [rb'github_pat_[A-Za-z0-9_]{20,}', rb'gh[pousr]_[A-Za-z0-9]{20,}',
                rb'sk-ant-[A-Za-z0-9_-]{20,}', rb'sk-[A-Za-z0-9_-]{32,}',
                rb'-----BEGIN (?:OPENSSH |RSA |EC )?PRIVATE KEY-----',
                rb'(?i)Bearer\s+[A-Za-z0-9_\-.]{24,}']
    rows, findings = [], []
    for file in sorted(set(allowed)):
        if not file.is_file() or file.is_symlink():
            continue
        data = file.read_bytes()
        relative = str(file.relative_to(stage))
        hit = any(re.search(pattern, data) for pattern in patterns)
        if hit:
            findings.append({'path': relative, 'reason': 'possible credential; content deliberately not printed'})
        rows.append({'path': relative, 'bytes': len(data), 'sha256': digest_file(file), 'sensitive_match': hit})
    output = run / 'reports'
    result = {'created_utc': utcnow(), 'scan_pass': not findings, 'findings': findings, 'files': rows,
              'limitations': 'Pattern scan, not a guarantee against all secrets; manual review remains required.',
              'excluded': ['cache/', 'snapshots/', 'task field arrays', 'credentials', 'complete sessions', '.git/'],
              'git_mutation_performed': False}
    atomic_json(output / 'github_upload_scan.json', result)
    (output / 'github_upload_allowlist.txt').write_text('\n'.join(row['path'] for row in rows if not row['sensitive_match']) + '\n')
    print(json.dumps({'scan_pass': not findings, 'file_count': len(rows), 'finding_count': len(findings),
                      'list': str(output / 'github_upload_allowlist.txt')}))
    return 0 if not findings else 1


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--run', required=True)
    args = parser.parse_args()
    raise SystemExit(main(args.run))
