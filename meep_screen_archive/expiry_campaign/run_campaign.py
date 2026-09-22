"""Detached finite queue; always package scientific raw evidence after exit."""
from pathlib import Path
import os
import subprocess
import sys
import time
from ti2d.queue import atomic_json

ROOT=Path(__file__).resolve().parent


def main():
    if sys.platform!='linux' or str(ROOT)!='/root/autodl-tmp/meep_sim/meep_screen/expiry_campaign':
        raise RuntimeError('SERVER_ONLY')
    os.chdir(ROOT)
    manifest=ROOT/'runs/expiry_20260922_H/manifest.json'
    code=1
    try:
        code=subprocess.call([sys.executable,'-B','-m','ti2d.queue','--manifest',str(manifest),'--phase','validation'])
    finally:
        export=subprocess.run([sys.executable,'-B',str(ROOT/'export_archive.py')],timeout=300).returncode
        # This archive never includes account credentials or unrelated legacy directories.
        out=ROOT/'exports/H_raw_evidence.tar.gz'
        temp=ROOT/'exports/H_raw_evidence.partial.tar.gz'
        raw=subprocess.run(['tar','-czf',str(temp),'-C',str(ROOT),'ti2d','runs','tests',
                            'create_campaign.py','export_archive.py','run_campaign.py','HANDOFF.md'],timeout=3600).returncode
        if raw==0:
            os.replace(temp,out)
            with (ROOT/'exports/H_raw_evidence.sha256').open('w') as f:
                subprocess.run(['sha256sum',str(out)],stdout=f,check=True)
        atomic_json(ROOT/'exports/final_status.json',dict(queue_exit_code=code,compact_export_exit_code=export,
            raw_archive_exit_code=raw,finished_epoch=time.time(),requires_off_server_backup=True,
            note='No automatic GitHub push from server; local archive branch update needs authenticated client.'))
    return code


if __name__=='__main__':
    raise SystemExit(main())
