import copy
import json
import math
import unittest

import numpy as np

from floquet_reconstruct import (FRACTIONS, analyze_case, closure_diagnostic,
                                maxwell_fields, parse_export, poynting_z)


def fixture(qx=0., qy=0., up=(.3+.1j,.2-.2j), down=(1.,.4j), extra=(), lambda_um=10.5, height_um=80.):
    lam=lambda_um*1e-6
    k0=2*math.pi/lam
    area=1e-9
    Z0=376.730313668
    height=height_um*1e-6
    cap=10.5e-6
    zref=height+.5*cap
    modes=[(0,0,qx*k0,qy*k0,np.array(up),np.array(down)),*extra]
    planes=[]
    for frac in FRACTIONS:
        z=height+frac*cap
        orders=[]
        for n1,n2,x,y,u,d in modes:
            kz2=k0*k0-x*x-y*y
            kz=complex(math.sqrt(kz2),0) if kz2>=0 else complex(0,-math.sqrt(-kz2))
            e=u*np.exp(-1j*kz*(z-zref))+d*np.exp(1j*kz*(z-zref))
            orders.append(dict(n1=n1,n2=n2,qx_per_m=x,qy_per_m=y,kz2_per_m2=kz2,Ex=[float(e[0].real),float(e[0].imag)],Ey=[float(e[1].real),float(e[1].imag)],padding11=bool(math.hypot(x,y)<=1.1*k0),padding12=True,prop=bool(kz2>0)))
        planes.append(dict(z_m=z,area_m2=area,orders=orders))
    return dict(case_id="synthetic",lambda_m=lam,Z0=Z0,nominal_Pin_W=1.,height_m=height,cap_m=cap,cap_fractions=list(FRACTIONS),planes=planes)


def lattice_fixture(theta_deg=15.,phi_deg=30.,lambda_um=10.5,height_um=80.):
    lam=lambda_um*1e-6
    k0=2*math.pi/lam
    pitch=80e-6
    bx=2*math.pi/pitch
    by=bx/math.sqrt(3)
    B=np.array([[bx,-by],[0,2*by]])
    theta,phi=math.radians(theta_deg),math.radians(phi_deg)
    kin=np.array([k0*math.sin(theta)*math.cos(phi),k0*math.sin(theta)*math.sin(phi)])
    area=math.sqrt(3)*pitch*pitch/2
    amp=math.sqrt(2*376.730313668*math.cos(theta)/area)
    down=np.array([amp*math.cos(phi),amp*math.sin(phi)])
    extra=[]
    for n1 in range(-25,26):
        for n2 in range(-25,26):
            if (n1,n2)==(0,0):
                continue
            q=kin+np.array([n1,n2])@B
            if np.linalg.norm(q)<=1.2*k0*(1+1e-10):
                extra.append((n1,n2,float(q[0]),float(q[1]),np.array([.2+.1j,.1]),np.array([0.,0.])))
    c=fixture(qx=kin[0]/k0,qy=kin[1]/k0,up=.4*down,down=down,extra=extra,lambda_um=lambda_um,height_um=height_um)
    c.update(reciprocal_basis_per_m=B.tolist(),kin_rad_per_m=[float(kin[0]),float(kin[1]),-k0*math.cos(theta)],expected_cell_area_m2=area)
    for p in c["planes"]:
        p["area_m2"]=area
    return c


class ReconstructionTests(unittest.TestCase):
    def assert_closed(self, result):
        for shell in result["fixed_nested_cutoffs"]:
            self.assertTrue(shell["field_prediction_checks_passed"])
            self.assertLess(max(shell["LOO_held_out_relative_L2_E_residual"]),1e-12)

    def test_normal_coherent_up_down(self):
        case=fixture()
        out=analyze_case(case)
        self.assert_closed(out)
        mode=out["mode_results"][0]
        self.assertTrue(np.allclose(np.array(mode["up_Ex_Ey_re_im"]),[[.3,.1],[.2,-.2]],atol=1e-12))
        self.assertTrue(np.allclose(np.array(mode["down_Ex_Ey_re_im"]),[[1.,0],[0,.4]],atol=1e-12))
        expected=case["planes"][0]["area_m2"]/(2*case["Z0"])*(abs(.3+.1j)**2+abs(.2-.2j)**2-1-.4**2)
        self.assertAlmostEqual(mode["planes"][0]["fitted_EH_net_W"],expected,places=22)

    def test_oblique_p_s_admittance_and_sign(self):
        for pol in ("p","s"):
            theta=math.radians(60)
            et=(1.,0.) if pol=="p" else (0.,1.)
            case=fixture(qx=math.sin(theta),up=(0.,0.),down=et)
            out=analyze_case(case)
            self.assert_closed(out)
            expect=case["planes"][0]["area_m2"]/(2*case["Z0"])*(1/math.cos(theta) if pol=="p" else math.cos(theta))
            self.assertAlmostEqual(out["fixed_nested_cutoffs"][0]["total_propagating_incoming_W"],expect,places=22)

    def test_rotated_oblique_p_s(self):
        phi=.4
        k0=2*math.pi/10.5e-6
        for e in ((math.cos(phi),math.sin(phi)),(-math.sin(phi),math.cos(phi))):
            out=analyze_case(fixture(qx=.5*math.cos(phi),qy=.5*math.sin(phi),up=(.2*e[0],.2*e[1]),down=e))
            self.assert_closed(out)
            shell=out["fixed_nested_cutoffs"][0]
            self.assertAlmostEqual(shell["propagating_outgoing_W"]/shell["total_propagating_incoming_W"],.04,places=12)

    def test_evanescent_cross_terms_preserved(self):
        k0=2*math.pi/10.5e-6
        extra=[(1,0,1.05*k0,0.,np.array([0.,1.]),np.array([0.,1j]))]
        out=analyze_case(fixture(extra=extra))
        self.assert_closed(out)
        ev=out["mode_results"][1]
        self.assertLess(ev["kz_re_im"][1],0)
        self.assertGreater(abs(ev["planes"][2]["up_down_cross_W"]),1e-15)
        self.assertAlmostEqual(ev["planes"][0]["fitted_EH_net_W"],ev["planes"][4]["fitted_EH_net_W"],places=22)
        self.assertAlmostEqual(ev["planes"][2]["up_self_W"],0.,places=22)

    def test_nonzero_incoming_not_forced_zero(self):
        k0=2*math.pi/10.5e-6
        extra=[(1,0,.3*k0,0.,np.array([.1,0]),np.array([.2,0]))]
        out=analyze_case(fixture(extra=extra))
        self.assertGreater(out["fixed_nested_cutoffs"][0]["nonzero_order_incoming_W_preserved"],0)

    def test_fixed_nested_cutoffs_not_cherry_picked(self):
        k0=2*math.pi/10.5e-6
        extra=[(1,0,1.15*k0,0.,np.array([0.,1.]),np.array([0.,1j]))]
        out=analyze_case(fixture(extra=extra))
        self.assertEqual([x["mode_count"] for x in out["fixed_nested_cutoffs"]],[1,2])
        self.assertNotEqual(out["fixed_nested_cutoffs"][0]["all_five_fitted_EH_net_W"],out["fixed_nested_cutoffs"][1]["all_five_fitted_EH_net_W"])

    def test_grazing_is_not_silently_dropped(self):
        out=analyze_case(fixture(qx=1.))
        self.assertEqual(len(out["mode_results"]),1)
        self.assertEqual(out["unresolved_orders_not_dropped"],[(0,0)])
        self.assertFalse(out["fixed_nested_cutoffs"][0]["complete_no_mode_silently_removed"])

    def test_bad_dispersion_rejected(self):
        c=fixture()
        c["planes"][2]["orders"][0]["kz2_per_m2"]*=.9
        with self.assertRaisesRegex(ValueError,"dispersion"):
            analyze_case(c)

    def test_inconsistent_order_inventory_rejected(self):
        c=fixture()
        c["planes"][2]["orders"]=[]
        with self.assertRaises(ValueError):
            analyze_case(c)

    def test_malformed_complex_and_incorrect_planes_rejected(self):
        for alteration in (lambda c:c["planes"][0]["orders"][0].update(Ex=[1]),lambda c:c.update(cap_fractions=[.1,.2,.5,.75,.9]),lambda c:c.update(nominal_Pin_W=.96)):
            c=fixture()
            alteration(c)
            with self.assertRaises(ValueError):
                analyze_case(c)

    def test_actual_held_out_error_detected(self):
        c=fixture()
        c["planes"][4]["orders"][0]["Ex"][0]+=.2
        out=analyze_case(c)
        self.assertFalse(out["fixed_nested_cutoffs"][0]["field_prediction_checks_passed"])
        self.assertGreater(max(out["fixed_nested_cutoffs"][0]["LOO_held_out_relative_L2_E_residual"]),.005)
        self.assertTrue(out["reconstructed_plane_invariance_is_algebraic_not_validation"])

    def test_no_acceptance_or_calibrating_loss(self):
        out=analyze_case(fixture())
        before=json.dumps(out,sort_keys=True)
        comp=closure_diagnostic(out,0.2)
        self.assertEqual(before,json.dumps(out,sort_keys=True))
        self.assertFalse(out["formal_qualified"])
        self.assertFalse(comp["loss_used_for_fit_or_normalization"])
        self.assertEqual(out["preserved_absolute_gate"],.005)
        self.assertEqual(out["preserved_relative_gate"],.02)

    def test_stream_and_single_case_parse(self):
        c=fixture()
        self.assertEqual(parse_export(json.dumps(c)),[c])
        self.assertEqual(parse_export(json.dumps({"cases":[c]})),[c])
        meta={k:v for k,v in c.items() if k!="planes"}
        lines=["COMSOL boot", "CODEX_FLOQUET_CASE_BEGIN_JSON "+json.dumps(meta)]
        lines.extend("CODEX_FLOQUET_PLANE_JSON "+json.dumps(dict(case_id=c["case_id"],**p)) for p in c["planes"])
        parsed=parse_export("\n".join(lines))
        self.assert_closed(analyze_case(parsed[0]))
        self.assert_closed(analyze_case(parse_export("\n".join(lines).replace("JSON ","JSON="))[0]))

    def test_independent_complete_lattice_and_1W(self):
        for theta,phi,prop,total in ((15.,30.,157,230),(60.,15.,156,225)):
            out=analyze_case(lattice_fixture(theta,phi))
            self.assert_closed(out)
            proof=out["independent_inventory_and_area_proof"]
            self.assertTrue(proof["verified"])
            self.assertEqual(proof["expected_total_orders"],total)
            self.assertEqual(proof["expected_propagating_orders"],prop)
            for s in out["fixed_nested_cutoffs"]:
                self.assertTrue(s["complete_no_mode_silently_removed"])
                self.assertTrue(s["specified_1W_incident_power_check_passed"])
                self.assertAlmostEqual(s["total_propagating_incoming_W"],1.,places=12)
            json.dumps(out,allow_nan=False)

    def test_same_missing_order_in_all_planes_rejected(self):
        c=lattice_fixture()
        for p in c["planes"]:
            p["orders"].pop()
        with self.assertRaisesRegex(ValueError,"incomplete independent"):
            analyze_case(c)

    def test_uniform_cut_area_and_shifted_bloch_rejected(self):
        c=lattice_fixture()
        for p in c["planes"]:
            p["area_m2"]*=.9
        with self.assertRaisesRegex(ValueError,"primitive-cell area"):
            analyze_case(c)
        c=lattice_fixture()
        for p in c["planes"]:
            p["orders"][1]["qx_per_m"]+=.2
            r=p["orders"][1]
            r["kz2_per_m2"]=(2*math.pi/c["lambda_m"])**2-r["qx_per_m"]**2-r["qy_per_m"]**2
        with self.assertRaisesRegex(ValueError,"fixed-kin reciprocal"):
            analyze_case(c)

    def test_missing_basis_never_claims_complete(self):
        out=analyze_case(fixture())
        self.assertFalse(out["independent_inventory_and_area_proof"]["verified"])
        self.assertFalse(out["fixed_nested_cutoffs"][0]["complete_no_mode_silently_removed"])

    def test_failed_markers_and_duplicate_cases(self):
        with self.assertRaises(ValueError):
            parse_export("no results")
        c=fixture()
        line="CODEX_FLOQUET_CASE_BEGIN_JSON "+json.dumps({k:v for k,v in c.items() if k!="planes"})
        with self.assertRaisesRegex(ValueError,"duplicate"):
            parse_export(line+"\n"+line)


if __name__=="__main__":
    unittest.main()

