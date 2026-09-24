"""Bounded Ordal-Ti two-dimensional blind-slot campaign."""

import os
from pathlib import Path

# Scope plotting and temporary caches to the authorized data-disk stage only.
# This does not change the user's home, environment installation or other apps.
_stage = Path(__file__).resolve().parents[1]
if str(_stage) == '/root/autodl-tmp/meep_sim/meep_screen/expiry_campaign':
    for _name, _subdir in (('MPLCONFIGDIR', 'cache/matplotlib'), ('TMPDIR', 'tmp')):
        _directory = _stage / _subdir
        _directory.mkdir(parents=True, exist_ok=True)
        os.environ[_name] = str(_directory)
