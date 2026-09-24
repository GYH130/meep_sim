"""Read-only, fixture-specific preliminary diagnostic. No Meep or file writes.

Not a replacement for the locked solver's acceptance, normalization or energy
audit. Applies only to the saved uniform-air p-polarized 30-degree L fixture.
"""
import argparse
import hashlib
import io
import json
from pathlib import Path
import tarfile

import numpy as np


def inspect(archive):
    expected = '0e96c13a36eb92e8b5fd9d01117c3a247f4493cc38a771d6ce3427fc2c6a3096'
    digest = hashlib.sha256()
    with archive.open('rb') as stream:
        for block in iter(lambda: stream.read(1048576), b''):
            digest.update(block)
    if digest.hexdigest() != expected:
        raise ValueError('Use the preserved 20260924T101823822823Z L raw evidence archive')
    rows = []
    with tarfile.open(archive, 'r:gz') as tf:
        for suffix in ('1e8', '1e6'):
            task = 'L_narrow_candidate_' + suffix
            prefix = 'boundary_screen_L/runs/boundary_20260924_L/tasks/' + task + '/attempt_1/'
            result = json.load(tf.extractfile(prefix + 'result.json'))
            case = result['case_config']
            if not (case['case'] == 'flat' and case['polarization'] == 'p'
                    and case['geometry']['period_um'] == 2
                    and abs(case['emission_theta_deg']) == 30
                    and case['wavelength_um'] == 10.5):
                raise ValueError('Not the fixed uniform L diagnostic fixture')
            with np.load(io.BytesIO(tf.extractfile(prefix + 'monitor_planes.npz').read()),
                         allow_pickle=False) as fields:
                ex, hz = fields['total_Ex'], fields['total_Hz']
                weights = fields['raw_refl_weights_canonical_uniform']
                if not (ex.shape == hz.shape == weights.shape == (64,)):
                    raise ValueError('Unexpected field dimensions; no silent cropping')
                q = np.cos(np.deg2rad(30.))
                up = (hz - ex / q) / 2
                down = (hz + ex / q) / 2
                power_up = float(np.sum(weights * np.abs(up)**2))
                power_down = float(np.sum(weights * np.abs(down)**2))
                if not np.isfinite(power_up + power_down) or power_down <= 0:
                    raise ValueError('Invalid directional amplitude denominator')
                ratio = power_up / power_down
            rows.append(dict(task=task, preliminary_total_field_ratio=ratio,
                analytical_R=result['analytical']['R'],
                preliminary_abs_difference=abs(ratio-result['analytical']['R']),
                official_R=result['metrics']['R'], original_status=result['status'],
                status_changed=False, source_npz=prefix+'monitor_planes.npz'))
    return dict(kind='PRELIMINARY_OFFLINE_DIAGNOSTIC_NOT_QUALIFICATION', rows=rows,
        assumptions='Uniform p wave in air at 30 degrees; same-plane up/down impedance factors cancel.',
        limitations='Continuum plane-wave impedance; must still audit Yee registration, background subtraction, actual incident power, independent absorption and full-cell geometry.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('raw_archive', type=Path)
    print(json.dumps(inspect(parser.parse_args().raw_archive), indent=2))
