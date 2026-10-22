"""Offline archive integrity checks; standard library only, no solver or network."""
import hashlib
import json
from pathlib import Path

root = Path(__file__).resolve().parents[1]
manifest = json.loads((root/'MANIFEST_SHA256.json').read_text())
for entry in manifest['files']:
    p = root/entry['path']
    raw = p.read_bytes()
    assert len(raw) == entry['bytes'], entry['path']
    assert hashlib.sha256(raw).hexdigest() == entry['sha256'], entry['path']
data = json.loads((root/'data/floorfill_3d_conditions.json').read_text())
rows = data['records']
assert len(rows) == 24 and sum(r['accepted'] for r in rows) == 5
keys = data['condition_key']
assert len({tuple(r[k] for k in keys) for r in rows}) == 24
for r in rows:
    assert abs(r['A_port'] + r['R_total'] - 1) < 1e-12
    assert abs(r['A_loss_total'] - r['A_film'] - r['A_Ti_H']) < 1e-12
    assert abs(r['port_loss_abs'] - abs(r['A_port']-r['A_loss_total'])) < 1e-12
    assert abs(r['flux_loss_abs'] - abs(r['A_flux']-r['A_loss_total'])) < 1e-12
    assert not r['formal_qualification'] and not r['mesh_convergence_established']
    assert r['solver_completed']
latest = json.loads((root/'data/latest_eight_quality.json').read_text())['cases']
assert len(latest) == 8 and sum(c['evaluation']['accepted'] for c in latest) == 1
assert all(c['native_exit_code'] == 0 for c in latest)
assert all(a['matches_parent'] for c in latest for a in c['artifact_checks'])
byid = {r['case_id']:r for r in rows}
for c in latest:
    r=byid[c['metrics']['case_id']]
    assert r['accepted'] == c['evaluation']['accepted']
    for k in ['A_port','A_loss_total','A_flux','A_film','A_Ti_H']:
        assert r[k] == c['metrics'][k]
assert len(list((root/'figures').glob('*.png'))) == 12
assert len(list((root/'report').glob('*.docx'))) == 1
assert len(list((root/'report').glob('*.pdf'))) == 1
print(json.dumps({'state':'VERIFIED_PUBLIC_RESEARCH_ARCHIVE','files':len(manifest['files']),
                  'unique_3d_conditions':24,'exploratory_accepted':5,'flux_warning':19,
                  'figures':12,'formal_qualification':False},ensure_ascii=False))
