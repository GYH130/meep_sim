"""最基本的烟雾测试：保证骨架可以 import，单位换算自洽。"""

import math

import numpy as np
import pytest

from src import geometry, io_utils, materials, postprocess, simulation  # noqa: F401


def test_wavelength_freq_roundtrip():
    for lam in [3.0, 8.0, 10.5, 13.0, 25.0]:
        f = materials.wavelength_um_to_meep_freq(lam)
        assert math.isclose(materials.meep_freq_to_wavelength_um(f), lam, rel_tol=1e-12)


def test_freq_range_for_band_ordering():
    f_min, f_max, f_c = materials.freq_range_for_band(8.0, 13.0)
    assert f_min < f_c < f_max
    assert math.isclose(f_max, 1.0 / 8.0)
    assert math.isclose(f_min, 1.0 / 13.0)


def test_wavelength_invalid():
    with pytest.raises(ValueError):
        materials.wavelength_um_to_meep_freq(0)
    with pytest.raises(ValueError):
        materials.freq_range_for_band(13.0, 8.0)


def test_estimate_resolution_positive():
    spec = geometry.GeometrySpec(
        kind="flat", period_x=0, period_y=0,
        feature_size=0.5, height=1.0, extra={},
    )
    res = geometry.estimate_resolution(spec, wavelength_min_um=8.0)
    assert res > 0
    assert isinstance(res, int)


def test_compute_RTA_energy_conservation():
    # 构造一个完美一致的人工 flux：R=0.3, T=0.2 → A=0.5
    ref = np.array([1.0, 1.0, 1.0])
    refl = -0.3 * ref
    trans = 0.2 * ref
    rta = postprocess.compute_RTA(refl, trans, ref)
    np.testing.assert_allclose(rta["R"], 0.3)
    np.testing.assert_allclose(rta["T"], 0.2)
    np.testing.assert_allclose(rta["A"], 0.5)
    report = postprocess.energy_conservation_check(rta)
    assert report["max_abs_dev"] < 1e-12


def test_compute_RTA_shape_mismatch():
    with pytest.raises(ValueError):
        postprocess.compute_RTA(np.zeros(3), np.zeros(4), np.ones(3))


def test_band_average_basic():
    wl = np.linspace(5.0, 15.0, 11)  # 5,6,...,15
    vals = np.ones_like(wl) * 0.8
    avg = postprocess.band_average(vals, wl, 8.0, 13.0)
    assert math.isclose(avg, 0.8)


def test_make_run_tag_format():
    tag = io_utils.make_run_tag("scripts/foo.py")
    assert tag.startswith("foo_")
    assert len(tag.split("_")[-1]) == 15  # YYYYmmdd-HHMMSS
