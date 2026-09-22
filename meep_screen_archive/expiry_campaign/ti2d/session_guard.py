"""Owned-session shutdown, copied from the already tested G diagnostic guard."""
import os
from pathlib import Path
import signal
import time


def _process_identity(pid, proc_root="/proc"):
    try:
        fields = (Path(proc_root) / str(pid) / "stat").read_text().rsplit(")", 1)[1].split()
    except (FileNotFoundError, ProcessLookupError):
        return None
    return {"pid": int(pid), "state": fields[0], "pgid": int(fields[2]),
            "session": int(fields[3]), "start_ticks": int(fields[19])}


def session_live_members(session_id, proc_root="/proc"):
    root = Path(proc_root)
    if not root.is_dir():
        raise RuntimeError("PROCESS_SESSION_AUDIT_UNAVAILABLE")
    members = []
    for path in root.iterdir():
        if not path.name.isdigit():
            continue
        record = _process_identity(int(path.name), root)
        if record and record["session"] == session_id and record["state"] not in ("Z", "X", "x"):
            members.append(record)
    return members


def signal_session_member(member, session_id, sig):
    current = _process_identity(member["pid"])
    if (current is None or current["session"] != session_id or
            current["start_ticks"] != member["start_ticks"] or current["state"] in ("Z", "X", "x")):
        return
    try:
        os.kill(current["pid"], sig)
    except ProcessLookupError:
        pass


def stop_owned_session(child, deadline, publish):
    started, errors, sent = time.monotonic(), set(), set()
    kill_at = min(started + 20., max(started, deadline - 5.))
    next_write = started
    while True:
        child.poll()
        try:
            members = session_live_members(child.pid)
        except Exception as exc:
            return {"stopped": False, "errors": [str(exc)], "live_members": None}
        if child.returncode is not None and not members:
            return {"stopped": True, "live_members": [], "errors": sorted(errors)}
        now = time.monotonic()
        if now >= deadline:
            return {"stopped": False, "live_members": members, "errors": sorted(errors)}
        sig = signal.SIGKILL if now >= kill_at else signal.SIGTERM
        for member in members:
            if time.monotonic() >= deadline:
                break
            identity = (member['pid'], member['start_ticks'], int(sig))
            if identity not in sent:
                try:
                    signal_session_member(member, child.pid, sig)
                    sent.add(identity)
                except Exception as exc:
                    errors.add(str(exc))
        if now >= next_write:
            publish()
            next_write = now + 5.
        time.sleep(min(.2, max(0., deadline-time.monotonic())))
