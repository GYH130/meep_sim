"""Allowlisted compact offline export; scans secrets; never uploads or runs Meep."""
from pathlib import Path
import hashlib
import json
import os
import re
import shutil
import tarfile
import tempfile
import time

ROOT = Path(__file__).resolve().parent
STAGE = ROOT.parent
PATTERNS = [rb'github_pat_[A-Za-z0-9_]{25,}', rb'gh[pousr]_[A-Za-z0-9]{25,}',
            rb'sk-(?:ant-|proj-)[A-Za-z0-9_-]{20,}',
            rb'-----BEGIN (?:OPENSSH |RSA |EC )?PRIVATE KEY-----',
            rb'https?://[^\s/:]+:[^\s/@]+@']


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def main():
    exports = ROOT/'exports'
    exports.mkdir(parents=True,exist_ok=True)
    target = exports/'compact_latest.tar.gz'
    with tempfile.TemporaryDirectory(prefix='archive_',dir=exports) as tmp:
        dest = Path(tmp)/'meep_screen_archive'
        dest.mkdir()
        inventory = []
        def add(path, relative=None):
            if not path.is_file() or path.is_symlink():
                return
            if path.stat().st_size > 8*1024*1024:
                raise ValueError('COMPACT_FILE_TOO_LARGE:' + str(path.relative_to(STAGE)))
            data=path.read_bytes()
            if any(re.search(pattern,data) for pattern in PATTERNS):
                raise ValueError('SECRET_SCAN_BLOCKED:' + str(path.relative_to(STAGE)))
            rel = relative or path.relative_to(STAGE)
            out=dest/rel
            out.parent.mkdir(parents=True,exist_ok=True)
            shutil.copyfile(path,out)
            inventory.append({'path':str(rel),'sha256':digest(out),'bytes':len(data)})
        for folder in (STAGE/'ti2d',STAGE/'assets',STAGE/'tests',ROOT/'ti2d',ROOT/'tests'):
            for path in sorted(folder.rglob('*')):
                if '__pycache__' not in path.parts and path.suffix in ('.py','.md','.json','.csv','.txt','.yml','.yaml'):
                    add(path)
        for path in sorted(ROOT.glob('*.py')): add(path)
        for path in sorted(STAGE.glob('*HANDOFF.md')): add(path)
        for path in sorted(ROOT.glob('*.md')): add(path)
        add(STAGE/'budget_ledger.json')
        sources=[STAGE/'runs/flat_r32_20260922_D', STAGE/'runs/baseline_r32_20260922_E',
                 STAGE/'repair_diagnostics/slot_stability_20260922_F',
                 STAGE/'repair_diagnostics/slot_stability_20260922_G_off', ROOT/'runs']
        verified, unqualified = [], []
        for source in sources:
            for path in sorted(source.rglob('result.json')):
                result=json.loads(path.read_text())
                if result.get('status')=='QUALIFIED':
                    evidence=result.get('evidence_sha256',{})
                    if not evidence:
                        raise ValueError('QUALIFIED_WITHOUT_EVIDENCE:'+str(path))
                    for name, expected in evidence.items():
                        p=path.parent/name
                        if p.is_symlink() or not p.resolve().is_relative_to(path.parent.resolve()) or digest(p)!=expected:
                            raise ValueError('EVIDENCE_HASH_FAILED:'+str(path))
                    verified.append({'path':str(path.relative_to(STAGE)),
                                     'status':'QUALIFIED_SINGLE_CASE_EVIDENCE_VERIFIED',
                                     'evidence_sha256':evidence})
                else:
                    unqualified.append({'path':str(path.relative_to(STAGE)),'status':result.get('status')})
            for path in sorted(source.rglob('*')):
                if path.name.startswith('.') or '__pycache__' in path.parts:
                    continue
                if path.suffix in ('.json','.md','.csv','.png','.svg') and not path.name.startswith('heartbeat'):
                    add(path)
        readme = '''# Ti 2D stage archive — incomplete study

Pure Ti, specified frozen Ordal Route A optical model; 10.5 µm, 2D only.
This is NOT a validated size–performance study, 3D-hole prediction or real high-temperature sample measurement.

## Evidence classification
- `VERIFICATION.json` lists individually qualified cases with checked raw-evidence SHA256.
- Other results remain failures, time-limited or diagnostic-only; no automatic upgrade.
- Old flat r32 with smoothing passed; old slanted r32 overflowed.
- Smoothing-off short test reached t≈98.83 before its wall limit; source ends at 525. It does NOT prove full stability.
- H is a 10-condition validation subset with smoothing off, serial 4 MPI ranks, fixed gates, no retries or size scan.
- A failed required condition stops expansion. Missing checks stay missing. Each case ≤6h; absolute stop 2026-09-23 08:52:44 UTC (16:52:44 China).

## Preservation and costs
This compact archive contains code, complete material inputs, configurations, small results/reports and raw-field fingerprints.
Large NPZ/HDF5 fields, caches, credentials, full conversation and raw logs are deliberately NOT committed.
Full pre-expiry raw evidence was separately downloaded as `preexpiry_evidence_20260922.tar.gz` (SHA256 17b1bc4eae4388bdcd1ccfe15496097ce20305102d0d3e5de1ab4eb8b1a24573).
An independent post-run raw bundle is generated on the server. Do not delete/release the instance before retrieving it.
`budget_ledger.json` measures union solver time and ranks; AI usage/billing and RMB costs are unavailable, not inferred. Background queue has no AI/API calls.

## Running / recovery
Remote stage `/root/autodl-tmp/meep_sim/meep_screen`, existing `meep_sim` Conda only.
See `expiry_campaign/HANDOFF.md`. Preserve absolute inputs or regenerate a NEW hashed manifest for another server; do not modify an already locked manifest.
This archive is a snapshot in a NEW GitHub branch. Main and the existing research branch are unchanged.
'''
        (dest/'README.md').write_text(readme)
        (dest/'VERIFICATION.json').write_text(json.dumps({'generated_epoch':time.time(),'scope':'individual cases only; full production gate remains incomplete','qualified':verified,'other_results':unqualified},indent=2)+'\n')
        (dest/'FILES.json').write_text(json.dumps(inventory,indent=2)+'\n')
        for path in dest.rglob('*'):
            if path.is_file() and any(re.search(p,path.read_bytes()) for p in PATTERNS):
                raise ValueError('FINAL_SECRET_SCAN_BLOCKED:'+str(path.relative_to(dest)))
        temporary=exports/('compact.tmp.'+str(os.getpid())+'.tar.gz')
        with tarfile.open(temporary,'w:gz') as tf:
            tf.add(dest,arcname=dest.name)
        os.replace(temporary,target)
        (exports/'compact_latest.sha256').write_text(digest(target)+'  compact_latest.tar.gz\n')
        print(json.dumps({'archive':str(target),'sha256':digest(target),'files':len(inventory),'qualified_cases':len(verified)}))


if __name__=='__main__':
    main()
