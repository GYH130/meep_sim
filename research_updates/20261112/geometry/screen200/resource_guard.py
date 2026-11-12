"""Resource checks for the resumed ideal-ray screen only.

User-approved admission: 72 GiB host RAM. Runtime reserve remains 64 GiB;
only dispatch concurrency is adapted. Owns explicitly retained Popen handles."""
from __future__ import annotations
import base64
import datetime as dt
import json
import math
import os
from pathlib import Path
import signal
import subprocess
import time

GB=1024**3
PS=Path(os.environ.get('MEEP_WINDOWS_POWERSHELL', 'powershell.exe'))
QUERY = r'''
$ErrorActionPreference='Stop'
$ProgressPreference='SilentlyContinue'
$os=Get-CimInstance Win32_OperatingSystem
$system=Get-CimInstance Win32_ComputerSystem
$volume=Get-Volume -DriveLetter D

# Win32_Processor sometimes omits CPU1's LoadPercentage on this host.
# Require performance counters covering all logical CPUs and each group's total.
$rows=@(Get-CimInstance Win32_PerfFormattedData_Counters_ProcessorInformation)
$logical=@($rows | Where-Object {$_.Name -match '^[0-9]+,[0-9]+$'})
$logicalNames=@($logical.Name | Sort-Object -Unique)
if($logicalNames.Count -ne $system.NumberOfLogicalProcessors -or $logical.Count -ne $logicalNames.Count -or $logical.Count -eq 0){
    throw 'Incomplete logical CPU counter coverage'
}
$groups=@($logical.Name | ForEach-Object {($_ -split ',')[0]} | Sort-Object -Unique)
$groupLoads=[ordered]@{}
foreach($group in $groups){
    $total=@($rows | Where-Object {$_.Name -eq ($group+',_Total')})
    if($total.Count -ne 1 -or $null -eq $total[0].PercentProcessorTime){
        throw ('Missing CPU group total '+$group)
    }
    $groupLoads[$group]=[double]$total[0].PercentProcessorTime
}
[pscustomobject]@{
    host_free_bytes=([double]$os.FreePhysicalMemory*1024)
    cpu_percent=([double](($groupLoads.Values | Measure-Object -Maximum).Maximum))
    cpu_group_loads=$groupLoads
    logical_cpu_counters=$logical.Count
    d_free_bytes=([double]$volume.SizeRemaining)
    d_health=$volume.HealthStatus.ToString()
}|ConvertTo-Json -Compress
'''
ENCODED_QUERY = base64.b64encode(QUERY.encode('utf-16le')).decode('ascii')

def rss(pid):
    try:
        for line in Path(f'/proc/{pid}/status').read_text().splitlines():
            if line.startswith('VmRSS:'):return int(line.split()[1])*1024
    except FileNotFoundError:return 0
    return 0

def owned_descendants(root_pids):
    parents={}
    for path in Path('/proc').glob('[0-9]*/stat'):
        try:
            rest=path.read_text().rsplit(')',1)[1].split()
            parents[int(path.parent.name)]=int(rest[1])
        except (FileNotFoundError,ProcessLookupError,PermissionError,ValueError,IndexError):continue
    selected=set(root_pids)
    while True:
        new={pid for pid,ppid in parents.items() if ppid in selected}
        if new<=selected:break
        selected.update(new)
    return selected

def admission_worker_limit(s):
    """Conservative dispatch cap, not a claim of CPU/RAM reservation.

    A lower cap drains existing workers naturally; it does not authorize killing
    them. The runner also caps this value by requested workers and CPU affinity.
    Runtime emergency stops still use the unchanged two-sample policy below.
    """
    if not (s.get('valid') is True and s['host_free_GiB']>=72
            and s['cpu_percent']<75 and s['D_free_GiB']>=100
            and s['D_health']=='Healthy' and s['wsl_available_GiB']>=8
            and s['owned_RSS_GiB']<112):
        return 0
    free,cpu=s['host_free_GiB'],s['cpu_percent']
    if free<80 or cpu>=65:
        return 2
    if free<100:
        return 4
    if free>=120 and cpu<55:
        return 16
    return 8

class Guard:
    def __init__(self):
        self.counters={'host_low':0,'owned_soft':0,'owned_hard':0,'wsl_low':0,'disk_bad':0}
    def sample(self,processes,controller_pid):
        now=dt.datetime.now(dt.timezone.utc).isoformat()
        # All passed processes were created by this runner. No process-name search.
        owned=[p for p in processes if p.poll() is None]
        try:
            q=subprocess.run([str(PS),'-NoLogo','-NoProfile','-NonInteractive','-EncodedCommand',ENCODED_QUERY],
                 capture_output=True,text=True,errors='replace',timeout=15,check=True)
            if len(q.stdout.encode())>32768 or len(q.stderr.encode())>32768:raise RuntimeError('Resource query exceeded 32KiB')
            host=json.loads(q.stdout.strip().lstrip('\ufeff'))
            for k in ('host_free_bytes','cpu_percent','d_free_bytes'):
                if not isinstance(host.get(k),(int,float)) or not math.isfinite(host[k]) or host[k]<0:raise RuntimeError('Missing resource '+k)
            if host['cpu_percent']>100:raise RuntimeError('Invalid CPU load')
            mi={k:int(v.split()[0])*1024 for k,v in (line.split(':',1) for line in Path('/proc/meminfo').read_text().splitlines())}
            tree=owned_descendants([p.pid for p in owned])
            s=dict(utc=now,valid=True,host_free_GiB=host['host_free_bytes']/GB,
              cpu_percent=host['cpu_percent'],cpu_group_loads=host['cpu_group_loads'],
              logical_cpu_counters=host['logical_cpu_counters'],D_free_GiB=host['d_free_bytes']/GB,D_health=host['d_health'],
              wsl_available_GiB=mi['MemAvailable']/GB,wsl_total_GiB=mi['MemTotal']/GB,
              owned_RSS_GiB=(rss(controller_pid)+sum(rss(pid) for pid in tree))/GB,
              owned_native_pids=[p.pid for p in owned],owned_process_tree_pids=sorted(tree))
            bad={'host_low':s['host_free_GiB']<64,'owned_soft':s['owned_RSS_GiB']>112,
                 'owned_hard':s['owned_RSS_GiB']>128,'wsl_low':s['wsl_available_GiB']<8,
                 'disk_bad':s['D_free_GiB']<100 or s['D_health']!='Healthy'}
            for k,v in bad.items():self.counters[k]=self.counters[k]+1 if v else 0
            s['consecutive_bad_samples']=dict(self.counters)
            s['stop_owned_reasons']=[k for k,v in self.counters.items() if v>=2]
            s['worker_limit']=admission_worker_limit(s)
            s['admit_new']=s['worker_limit']>0
            s['admission_host_min_GiB']=72
            s['runtime_host_min_GiB']=64
            return s
        except Exception as e:
            # A missing sample cannot authorize another child. It also breaks
            # consecutiveness, rather than fabricating a bad numeric reading.
            self.counters={k:0 for k in self.counters}
            details={}
            if isinstance(e,(subprocess.CalledProcessError,subprocess.TimeoutExpired)):
                for key in ('stdout','stderr'):
                    value=getattr(e,key,None)
                    if isinstance(value,bytes):value=value.decode('utf-8',errors='replace')
                    if value is not None:details[key]=str(value)[-4096:]
                details['returncode']=getattr(e,'returncode',None)
            return dict(utc=now,valid=False,error=repr(e),query_error_details=details,admit_new=False,worker_limit=0,admission_host_min_GiB=72,runtime_host_min_GiB=64,stop_owned_reasons=[])

def stop_owned(processes):
    """Stop only still-running children represented by retained Popen objects."""
    terminated=[]
    for p in processes:
        if p.poll() is None:
            # Each child was created with start_new_session=True. This process
            # group contains only its descendants, including compiler helpers.
            try:os.killpg(p.pid,signal.SIGTERM)
            except ProcessLookupError:continue
            terminated.append(p)
    deadline=time.monotonic()+5
    for p in terminated:
        try:p.wait(timeout=max(.01,deadline-time.monotonic()))
        except subprocess.TimeoutExpired:
            if p.poll() is None:
                try:os.killpg(p.pid,signal.SIGKILL)
                except ProcessLookupError:pass
            p.wait()
    return [p.pid for p in terminated]
