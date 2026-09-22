"""Server-only zero-step API/coordinate check; not a numerical qualification."""
import json
import meep as mp
import numpy as np

mp.verbosity(0)
s = mp.Simulation(cell_size=mp.Vector3(4, 8), resolution=4, Courant=.25,
                  k_point=mp.Vector3(.05, 0),
                  boundary_layers=[mp.PML(1, direction=mp.Y)],
                  sources=[mp.Source(mp.GaussianSource(.1, fwidth=.02), mp.Hz, center=mp.Vector3())])
vol = mp.Volume(center=mp.Vector3(), size=mp.Vector3(4, 4))
d = s.add_dft_fields([mp.Ex, mp.Dx, mp.Ey, mp.Dy], [.1], where=vol, yee_grid=True)
line = mp.Volume(center=mp.Vector3(0, 1), size=mp.Vector3(4, 0))
dl = s.add_dft_fields([mp.Ex, mp.Hz], [.1], where=line, yee_grid=False)
s.init_sim()
s.run(until=0)
out = {}
for c in (mp.Ex, mp.Dx, mp.Ey, mp.Dy):
    a = s.get_dft_array(d, c, 0)
    dims, lo, hi = s.get_array_slice_dimensions(c, vol=vol)
    out[mp.component_name(c)] = {'shape': a.shape, 'dims': dims.tolist(),
                                'lo': [lo.x, lo.y, lo.z], 'hi': [hi.x, hi.y, hi.z]}
for c in (mp.Ex, mp.Hz):
    a = s.get_dft_array(dl, c, 0)
    x, y, z, w = s.get_array_metadata(dft_cell=dl)
    out['line_'+mp.component_name(c)] = {'shape': a.shape, 'x': x.tolist(), 'y': y.tolist(), 'w_shape': w.shape, 'w_sum': float(w.sum())}
if mp.am_master():
    print(json.dumps(out))
s.reset_meep()
