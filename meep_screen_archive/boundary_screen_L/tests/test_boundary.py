"""Offline boundary selection/cache isolation; these tests never run FDTD."""
import ast
from copy import deepcopy
from pathlib import Path
from types import SimpleNamespace
import unittest

from ti2d.boundary import DEFAULT_PROFILE, boundary_spec, build_boundary_layers
from ti2d.common import digest_object
from ti2d.solver import layout, result_checks


def case(kind='absorber'):
    return {'boundary_kind': kind, 'boundary_profile': dict(DEFAULT_PROFILE), 'pml_um': 12.,
            'air_above_um': 12., 'air_below_um': 12.,
            'geometry': {'period_um': 50., 'thickness_um': 60.}}


class BoundaryTests(unittest.TestCase):
    def test_queue_reserves_evidence_time_even_when_deadline_shortens_case(self):
        source = (Path(__file__).resolve().parents[1] / 'ti2d/queue.py').read_text()
        tree = ast.parse(source)
        values = [k.value for n in ast.walk(tree) if isinstance(n, ast.Call)
                  for k in n.keywords if k.arg == 'TI2D_WALL_TIMEOUT_SECONDS']
        self.assertEqual(len(values), 1)
        expression = compile(ast.Expression(values[0]), '<timeout-expression>', 'eval')
        self.assertEqual(eval(expression, {'effective_timeout': 600., 'max': max, 'str': str}), '480.0')
        self.assertEqual(eval(expression, {'effective_timeout': 60., 'max': max, 'str': str}), '1.0')

    def test_explicit_absorber_factory_matches_installed_api(self):
        mp = SimpleNamespace(Y=1, ALL=-1,
                             Absorber=lambda **kw: ('absorber', kw),
                             PML=lambda **kw: ('pml', kw))
        name, kwargs = build_boundary_layers(mp, case())[0]
        self.assertEqual(name, 'absorber')
        self.assertEqual({k: v for k, v in kwargs.items() if k != 'pml_profile'},
                         {'thickness': 12., 'direction': 1, 'side': -1,
                          'R_asymptotic': 1e-15, 'mean_stretch': 1.})
        self.assertEqual([kwargs['pml_profile'](v) for v in (0., .5, 1.)], [0., .25, 1.])
        self.assertEqual(build_boundary_layers(mp, case('pml'))[0][0], 'pml')

    def test_legacy_pml_is_not_changed_to_absorber(self):
        config = {'pml_um': 12.}
        self.assertEqual(boundary_spec(config), boundary_spec(case('pml')))
        self.assertEqual(config, {'pml_um': 12.})

    def test_absorber_requires_explicit_profile_and_rejects_unknown_parameters(self):
        for changes in ({'boundary_kind': 'unknown'}, {'boundary_profile': None},
                        {'boundary_profile': {'name': 'quadratic'}},
                        {'boundary_profile': {**DEFAULT_PROFILE, 'name': 'linear'}},
                        {'boundary_profile': {**DEFAULT_PROFILE, 'mean_stretch': 2.}},
                        {'boundary_profile': {**DEFAULT_PROFILE, 'R_asymptotic': 1.}}):
            with self.subTest(changes=changes), self.assertRaises(ValueError):
                boundary_spec({**case(), **changes})
        for value in (0, -1, True, '12', float('nan'), float('inf')):
            with self.subTest(value=value), self.assertRaises(ValueError):
                boundary_spec({**case(), 'pml_um': value})

    def test_cache_boundary_identity_separates_kind_thickness_and_profile(self):
        base = boundary_spec(case())
        variations = [boundary_spec(case('pml')),
                      boundary_spec({**case(), 'pml_um': 6.}),
                      boundary_spec({**case(), 'boundary_profile':
                                     {**DEFAULT_PROFILE, 'R_asymptotic': 1e-12}})]
        for other in variations:
            self.assertNotEqual(digest_object({'boundary': base}), digest_object({'boundary': other}))
        source = (Path(__file__).resolve().parents[1] / 'ti2d/solver.py').read_text()
        tree = ast.parse(source)
        updates = [n for n in ast.walk(tree) if isinstance(n, ast.Call)
                   and isinstance(n.func, ast.Attribute) and n.func.attr == 'update'
                   and isinstance(n.func.value, ast.Name) and n.func.value.id == 'reference_identity']
        self.assertEqual(len(updates), 1)
        self.assertTrue(any(k.arg == 'boundary' for k in updates[0].keywords))
        self.assertIn('boundary=boundary', source)

    def test_boundary_type_does_not_change_physical_layout(self):
        a, p = case(), case('pml')
        self.assertEqual(layout(a), layout(p))
        self.assertEqual(layout(a)['height'], 108.)
        self.assertEqual(layout(a)['surface'], 30.)
        self.assertEqual(layout(a)['non_pml_height'], 84.)

    def test_accuracy_and_stop_thresholds_remain_unchanged(self):
        config = {**case(), 'case': 'flat'}
        metrics = {'R': .9, 'T': 0., 'A_flux': .1, 'A_vol': .1,
                   'order_R': .9, 'order_T': 0., 'reference_pseudo_reflection': 0.}
        analytic = {'R': .9, 'T': 0., 'A': .1}
        self.assertTrue(all(result_checks(config, metrics, True, analytic).values()))
        self.assertFalse(result_checks(config, metrics, True, analytic,
                                       'TIME_LIMIT_UNQUALIFIED')['stopped'])
        changed = deepcopy(metrics)
        changed['A_vol'] = .110001
        self.assertFalse(result_checks(config, changed, True, analytic)['independent_absorption'])
        self.assertFalse(result_checks(config, metrics, True,
                                       {**analytic, 'R': .902501})['fresnel_R'])


if __name__ == '__main__':
    unittest.main()
