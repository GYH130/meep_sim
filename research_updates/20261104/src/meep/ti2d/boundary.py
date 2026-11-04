"""Explicit, cache-keyed outer boundary; no Meep import during offline checks.

Meep 1.33.0 Absorber is a scalar electric and magnetic conductivity layer,
not a perfectly matched layer. Its API matches PML. The supported quadratic
profile and numerical defaults are made explicit to prevent cache mixing.
Physical layout retains the historical ``pml_um`` thickness key.
"""
from __future__ import annotations

import math


DEFAULT_PROFILE = {'name': 'quadratic', 'R_asymptotic': 1e-15, 'mean_stretch': 1.0}


def _positive_real(value, name):
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError('BOUNDARY_' + name + '_INVALID')
    if not math.isfinite(value) or value <= 0:
        raise ValueError('BOUNDARY_' + name + '_INVALID')
    return float(value)


def boundary_spec(config):
    """Return the exact constructor identity used for both model and reference.

    Missing kind means legacy PML. Absorber must be opted into explicitly and
    must carry its full profile config; no new absorbing layer is implicit.
    """
    kind = config.get('boundary_kind', 'pml')
    if kind not in ('pml', 'absorber'):
        raise ValueError('BOUNDARY_KIND_UNSUPPORTED')
    profile = config.get('boundary_profile')
    if profile is None:
        if kind == 'absorber':
            raise ValueError('ABSORBER_REQUIRES_EXPLICIT_BOUNDARY_PROFILE')
        profile = dict(DEFAULT_PROFILE)
    if not isinstance(profile, dict) or set(profile) != set(DEFAULT_PROFILE):
        raise ValueError('BOUNDARY_PROFILE_FIELDS_INVALID')
    if profile['name'] != 'quadratic':
        raise ValueError('BOUNDARY_PROFILE_UNSUPPORTED')
    reflectivity = _positive_real(profile['R_asymptotic'], 'R_ASYMPTOTIC')
    stretch = _positive_real(profile['mean_stretch'], 'MEAN_STRETCH')
    if reflectivity >= 1:
        raise ValueError('BOUNDARY_R_ASYMPTOTIC_INVALID')
    if kind == 'absorber' and stretch != 1.0:
        raise ValueError('ABSORBER_MEAN_STRETCH_MUST_EQUAL_ONE')
    return {'kind': kind, 'thickness_um': _positive_real(config['pml_um'], 'THICKNESS'),
            'direction': 'Y', 'side': 'ALL',
            'profile': {'name': 'quadratic', 'R_asymptotic': reflectivity,
                        'mean_stretch': stretch}}


def quadratic_profile(u):
    return u * u


def build_boundary_layers(mp, config):
    spec = boundary_spec(config)
    constructor = mp.Absorber if spec['kind'] == 'absorber' else mp.PML
    return [constructor(thickness=spec['thickness_um'], direction=mp.Y, side=mp.ALL,
                        R_asymptotic=spec['profile']['R_asymptotic'],
                        mean_stretch=spec['profile']['mean_stretch'],
                        pml_profile=quadratic_profile)]
