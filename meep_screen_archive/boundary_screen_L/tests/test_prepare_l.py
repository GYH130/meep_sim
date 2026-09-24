"""Pure input checks only; no solver launch or Meep import."""
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from prepare_l import (build_cases, control_proof, ROOT, CONTROL_ID, CANDIDATE_IDS,
                       FULL_IDS, K_INPUT_POWER, K_PSEUDO_REFLECTION, SOLVER_DEADLINE,
                       RENTAL_EXPIRY)
from ti2d.geometry import validate_geometry
from ti2d.solver import layout


def base_case():
    return dict(id='K_flat32p', case='flat', material_path='old.json', implementation_sha256='old',
                resolution=32, courant=.25, wavelength_um=10.5, polarization='p',
                emission_theta_deg=30, source_fwidth_fraction=.2, source_cutoff=5.,
                pml_um=12., air_above_um=12., air_below_um=12., sample_interval_um=2.,
                stop_window_um=20., stop_consecutive_windows=3, intensity_decay=1e-4,
                flux_drift_tolerance=.001, mpi_ranks=4, eps_averaging=False,
                geometry=dict(period_um=50., surface_fill=.692820323, axis_length_um=40.,
                              tilt_deg=30, thickness_um=60., min_wall_um=1., residual_um=10.))


class LPreparationTests(unittest.TestCase):
    def test_fixed_list_has_control_two_candidates_and_two_conditional_full_cases(self):
        cases = build_cases(base_case(), 'new.json', 'same_core')
        self.assertEqual([v[0] for v in cases], [CONTROL_ID] + CANDIDATE_IDS + FULL_IDS)
        self.assertEqual([v[2] for v in cases], [1200., 1200., 1200., 14400., 14400.])
        self.assertTrue(all(v[1]['case'] == 'flat' for v in cases))
        self.assertEqual([v[1]['boundary_profile']['R_asymptotic'] for v in cases],
                         [1e-15, 1e-8, 1e-6, 1e-8, 1e-6])

    def test_no_scientific_gate_source_material_or_y_layout_is_relaxed(self):
        base = base_case()
        for _, case, cap in build_cases(base, 'new.json', 'same_core'):
            for field in ('resolution','courant','wavelength_um','polarization','emission_theta_deg',
                          'source_fwidth_fraction','source_cutoff','sample_interval_um','stop_window_um',
                          'stop_consecutive_windows','intensity_decay','flux_drift_tolerance','mpi_ranks',
                          'eps_averaging'):
                self.assertEqual(case[field], base[field])
            self.assertTrue(validate_geometry(case['geometry'])['valid'])
            self.assertEqual({k:v for k,v in layout(case).items() if k!='period'},
                             {k:v for k,v in layout(base).items() if k!='period'})
            self.assertEqual(case['max_solver_seconds'], cap-120.)

    def test_placeholder_has_no_claim_to_be_a_size_design(self):
        cases = build_cases(base_case(), 'new.json', 'same_core')
        self.assertEqual([c['geometry']['period_um'] for _,c,_ in cases], [2.,2.,2.,50.,50.])
        self.assertTrue(all('not a scientific period' in c['diagnostic_scope'] for _,c,_ in cases[:3]))
        self.assertEqual(base_case()['geometry']['period_um'], 50.)

    def test_unchanged_core_bytes_match_K(self):
        origin = ROOT.parent / 'boundary_absorber_K/ti2d'
        for old in origin.glob('*.py'):
            self.assertEqual(old.read_bytes(), (ROOT/'ti2d'/old.name).read_bytes(), old.name)

    def test_control_proof_never_calls_expected_failure_qualified(self):
        checks = {k:True for k in ('stop','downward_flux','direction','plane_closure','empty_transmission')}
        checks['empty_pseudo_reflection'] = False
        result = dict(status='NUMERICAL_FAILED',failure_stage='reference',reference=dict(
            status='NUMERICAL_FAILED', checks=checks, identity={}, input_power=K_INPUT_POWER,
            pseudo_reflection=K_PSEUDO_REFLECTION, front_projection={},back_projection={}))
        with TemporaryDirectory() as d:
            path=Path(d)/'result.json'; path.write_text('{}')
            proof=control_proof(result,base_case(),path)
            self.assertEqual(proof['input_power_per_width'], K_INPUT_POWER/50)
            self.assertFalse(proof['control_comparison']['diagnostic_reproduction_is_qualification'])
            result['reference']['checks']['stop'] = False
            with self.assertRaisesRegex(ValueError, 'OTHER_CHECK_FAILED'):
                control_proof(result,base_case(),path)

    def test_deadline_reserves_83_minutes(self):
        self.assertEqual(RENTAL_EXPIRY-SOLVER_DEADLINE, 83*60)


if __name__ == '__main__':
    unittest.main()
