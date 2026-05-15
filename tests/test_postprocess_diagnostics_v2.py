import cmath
import math

import numpy as np
import pytest

from src import postprocess
from src.materials import get_ti_medium
from src.simulation import run_periodic_2d_metal_single_wavelength


def test_wavelength_integrated_average_uses_trapezoid_not_simple_mean():
    wavelengths = np.array([13.0, 8.0, 9.0])
    values = np.array([1.0, 0.0, 0.0])
    simple = float(np.mean(values))
    integrated = postprocess.wavelength_integrated_average(values, wavelengths, 8.0, 13.0)
    assert not math.isclose(integrated, simple)
    assert math.isclose(integrated, 0.4)
    assert postprocess.wavelength_integrated_average.last_status == "valid"


def test_downward_incidence_signed_transmission_definition():
    rta = postprocess.compute_RTA_downward_incidence(
        reflection_flux_raw=np.array([8.0]),
        transmission_flux_raw=np.array([0.42]),
        input_flux_raw=np.array([-10.0]),
    )
    np.testing.assert_allclose(rta["R"], [0.8])
    np.testing.assert_allclose(rta["raw_transmittance"], [0.042])
    np.testing.assert_allclose(rta["T"], [-0.042])
    np.testing.assert_allclose(rta["A"], [0.242])


def test_low_source_snr_is_flagged():
    quality = postprocess.broadband_source_quality_mask(
        np.array([1.0, 1e-4, 0.5]),
        relative_threshold=1e-3,
    )
    assert quality["source_quality_flag"].tolist() == [
        "VALID",
        "LOW_SOURCE_SNR",
        "VALID",
    ]
    assert quality["valid_for_quantitative_metric"].tolist() == [True, False, True]


def test_negative_transmission_failure_is_flagged():
    check = postprocess.opaque_substrate_transmission_check(np.array([-0.042]))
    assert check["transmission_quality_flag"].tolist() == ["FAIL"]
    assert check["valid_for_opaque_substrate_metric"].tolist() == [False]


def test_single_wavelength_flat_ti_10um_matches_fresnel():
    mp = pytest.importorskip("meep")
    ti = get_ti_medium(lambda_min_um=10.0, lambda_max_um=10.0)

    def factory(y_surface_um: float, substrate_thickness_um: float) -> list:
        return [
            mp.Block(
                material=ti,
                center=mp.Vector3(0, y_surface_um - 0.5 * substrate_thickness_um, 0),
                size=mp.Vector3(10.0, substrate_thickness_um, mp.inf),
            )
        ]

    result = run_periodic_2d_metal_single_wavelength(
        geometry_factory=factory,
        period_um=10.0,
        wavelength_um=10.0,
        resolution=32,
        pml_thickness_um=2.0,
        substrate_thickness_um=4.0,
        air_buffer_um=4.0,
        decay_db=20.0,
        source_component="Ez",
    )
    eps = ti.epsilon(0.1)[0][0]
    n_ti = cmath.sqrt(eps)
    r_fresnel = abs((1.0 - n_ti) / (1.0 + n_ti)) ** 2
    a_fresnel = 1.0 - r_fresnel
    assert abs(result["A"] - 0.0718) < 0.005
    assert abs(result["A"] - a_fresnel) < 0.01
    assert abs(result["T"]) < 1e-3
    assert result["source_mode"] == "single_wavelength"
    assert result["solver_version"] == "diagnostics_v2_signed_flux_single_wavelength"
