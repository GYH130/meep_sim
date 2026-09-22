"""One-shot completion-to-archive pipeline, no FDTD and no AI/API calls.

Run only for this campaign. Requires this Mac awake, SSH access, and its already
authorized GitHub SSH identity. No credential extraction, new keys or main writes.
"""
from datetime import datetime, timezone
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import tarfile
import tempfile
import time

BASE=Path(__file__).resolve().parent
REPO=BASE.parents[1]/'github_archive_20260922'
REMOTE='/root/autodl-tmp/meep_sim/meep_screen/expiry_campaign'
HOST='root@connect.nmb2.seetacloud.com'
SSH=['ssh','-o','BatchMode=yes','-o','ConnectTimeout=12','-o','ServerAliveInterval=30','-o','ServerAliveCountMax=3','-p','23417',HOST]
BRANCH='codex/archive-20260922-ti2d'
LOCAL_CUTOFF=datetime(2026,9,23,12,22,44,tzinfo=timezone.utc).timestamp()
READY_CUTOFF=datetime(2026,9,23,9,52,44,tzinfo=timezone.utc).timestamp()
SECRETS=[rb'github_pat_[A-Za-z0-9_]{25,}',rb'gh[pousr]_[A-Za-z0-9]{25,}',
         rb'sk-(?:ant-|proj-)[A-Za-z0-9_-]{20,}',rb'-----BEGIN (?:OPENSSH |RSA |EC )?PRIVATE KEY-----',
         rb'https?://[^\s/:]+:[^\s/@]+@']


def status(state,**kwargs):
    obj=dict(state=state,updated_utc=datetime.now(timezone.utc).isoformat(),**kwargs)
    tmp=BASE/'collector_status.tmp.json'
    tmp.write_text(json.dumps(obj,ensure_ascii=False,indent=2)+'\n')
    os.replace(tmp,BASE/'collector_status.json')
    print(json.dumps(obj,ensure_ascii=False),flush=True)


def run(args,timeout=90,check=True):
    left=LOCAL_CUTOFF-time.time()
    if left<=0: raise RuntimeError('LOCAL_RENTAL_SAFETY_CUTOFF')
    return subprocess.run(args,capture_output=True,text=True,check=check,
                          timeout=min(timeout,left),env=dict(os.environ,GIT_TERMINAL_PROMPT='0'))


def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''): h.update(b)
    return h.hexdigest()


def fetch(name,timeout=1800):
    path=BASE/name
    run(['rsync','-a','--partial','--timeout=120','-e',
         'ssh -p 23417 -o BatchMode=yes -o ConnectTimeout=12',
         HOST+':'+REMOTE+'/exports/'+name,str(path)],timeout=timeout)
    return path


def install_compact(path):
    if run(['git','-C',str(REPO),'status','--porcelain']).stdout.strip():
        raise RuntimeError('ARCHIVE_CHECKOUT_HAS_UNEXPECTED_LOCAL_EDITS')
    if run(['git','-C',str(REPO),'branch','--show-current']).stdout.strip()!=BRANCH:
        raise RuntimeError('ARCHIVE_BRANCH_CHANGED')
    with tarfile.open(path,'r:gz') as tf:
        members=tf.getmembers()
        for m in members:
            p=Path(m.name)
            if p.is_absolute() or '..' in p.parts or not p.parts or p.parts[0]!='meep_screen_archive':
                raise ValueError('UNSAFE_ARCHIVE_PATH')
            if not (m.isdir() or m.isfile()) or m.size>8*1024*1024:
                raise ValueError('UNSAFE_ARCHIVE_MEMBER')
        for m in members:
            if m.isfile() and any(re.search(p,tf.extractfile(m).read()) for p in SECRETS):
                raise ValueError('SECRET_SCAN_BLOCKED:'+m.name)
        # Archive paths have all been validated; write only this dedicated snapshot checkout.
        tf.extractall(REPO,members=members)
    inventory=json.loads((REPO/'meep_screen_archive/FILES.json').read_text())
    for item in inventory:
        p=(REPO/'meep_screen_archive'/item['path']).resolve()
        if not p.is_relative_to((REPO/'meep_screen_archive').resolve()) or sha(p)!=item['sha256']:
            raise ValueError('INVENTORY_HASH_MISMATCH')
    run(['git','-C',str(REPO),'add','--','meep_screen_archive'])
    run(['git','-C',str(REPO),'diff','--cached','--check'])
    if run(['git','-C',str(REPO),'diff','--cached','--quiet'],check=False).returncode:
        run(['git','-C',str(REPO),'commit','-m','Archive final bounded H evidence and qualification status'])
    run(['git','-C',str(REPO),'push','origin','HEAD:refs/heads/'+BRANCH],timeout=300)
    local=run(['git','-C',str(REPO),'rev-parse','HEAD']).stdout.strip()
    remote=run(['git','-C',str(REPO),'ls-remote','origin','refs/heads/'+BRANCH]).stdout.split()[0]
    if local!=remote: raise RuntimeError('REMOTE_COMMIT_NOT_VERIFIED')
    return local


def main():
    with (BASE/'collector.lock').open('a+') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        status('WAITING_FOR_SERVER_EXPORT',pid=os.getpid(),note='One-shot backup pipeline, no AI calls; Mac must stay awake and online')
        ready=None
        # Deterministic file-completion wait, not a model/tool polling conversation.
        while time.time()<READY_CUTOFF:
            try:
                r=run(SSH+['test -f '+REMOTE+'/exports/final_status.json && cat '+REMOTE+'/exports/final_status.json'],timeout=30,check=False)
                if r.returncode==0:
                    ready=json.loads(r.stdout)
                    break
            except (subprocess.SubprocessError,ValueError,OSError):
                pass
            time.sleep(min(600,max(0,READY_CUTOFF-time.time())))
        if ready is None:
            raise RuntimeError('NO_FINAL_SERVER_EXPORT_BEFORE_RESERVED_BACKUP_WINDOW')
        status('FETCHING_COMPACT',server_final_status=ready)
        checksum=fetch('compact_latest.sha256',90).read_text().split()[0]
        compact=fetch('compact_latest.tar.gz',900)
        if sha(compact)!=checksum: raise RuntimeError('COMPACT_TRANSFER_HASH_MISMATCH')
        commit=install_compact(compact)
        status('GITHUB_UPDATED_FETCHING_RAW',commit=commit)
        if ready.get('raw_archive_exit_code')!=0:
            raise RuntimeError('SERVER_RAW_EXPORT_FAILED_COMPACT_ALREADY_SAVED')
        expected=fetch('H_raw_evidence.sha256',90).read_text().split()[0]
        raw=fetch('H_raw_evidence.tar.gz',7200)
        if sha(raw)!=expected: raise RuntimeError('RAW_TRANSFER_HASH_MISMATCH')
        status('COMPLETE',commit=commit,raw_path=str(raw),raw_sha256=expected,
               branch_url='https://github.com/GYH130/meep_sim/tree/'+BRANCH,
               note='Compact GitHub snapshot and full raw H evidence verified; scientific qualification unchanged')


if __name__=='__main__':
    try:
        main()
    except Exception as exc:
        # Do not log credential-bearing environment or subprocess command output.
        status('BACKUP_NEEDS_ATTENTION',error_type=type(exc).__name__,message=str(exc)[:500])
        raise SystemExit(1)
