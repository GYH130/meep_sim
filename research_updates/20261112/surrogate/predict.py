"""Predict one scalar directional absorptivity. No FEM/server calls.

Example: python predict.py --height 90 --film 2.5 --theta 30 --phi 30 --pol s
New b or t beyond the training range requires --allow-extrapolation and remains diagnostic.
"""
from pathlib import Path
import argparse,json,sys,subprocess,time
import numpy as np, joblib
from scipy.spatial import Delaunay
from scipy.optimize import linprog
from surrogate import predict,features
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE/'ray_engine'))
import raylib

def ray_predict(row,rays=131072):
    args=raylib.build_args(row['height_um'],row['film_um'],row['wavelength_um'],
        row['theta_deg'],row['phi_deg'],rays,20261110,d_bottom=row['bottom_width_um'])
    raw=json.loads(subprocess.run([str(HERE/'ray_engine/ray_material')]+args,
                 check=True,capture_output=True,text=True,timeout=30).stdout)
    channel=raylib.channel(raw,row['polarization'])
    return channel

def evaluate(row,allow_extrapolation=False,method='gp',pool='primary_exploratory',rays=131072):
    start=time.monotonic()
    d=json.loads((HERE/'dataset.json').read_text());ids=set(d['pools'][pool])
    train=[r for r in d['rows'] if r['case_id'] in ids]
    tests={'bottom_width_seen':any(abs(r['bottom_width_um']-row['bottom_width_um'])<1e-9 for r in train),
      'wavelength_seen':any(abs(r['wavelength_um']-row['wavelength_um'])<1e-9 for r in train),
      'theta_seen':any(abs(r['theta_deg']-row['theta_deg'])<1e-9 for r in train),
      'phi_seen_or_normal':row['theta_deg']==0 or any(abs(r['phi_deg']-row['phi_deg'])<1e-9 for r in train)}
    x=np.unique(np.array([[r['height_um'],r['film_um']] for r in train]),axis=0)
    tests['inside_H_t_training_hull']=bool(Delaunay(x).find_simplex([row['height_um'],row['film_um']],tol=1e-9)>=0)
    # Marginal ranges alone miss combinations absent from the actual design table.
    x_joint=features(train);q_joint=features([row])[0]
    scale=np.where(np.ptp(x_joint,axis=0)>1e-10,np.ptp(x_joint,axis=0),1.)
    x_joint=(x_joint-x_joint[0])/scale;q_joint=(q_joint-features(train)[0])/scale
    joint=linprog(np.zeros(len(train)),A_eq=np.vstack([x_joint.T,np.ones(len(train))]),
        b_eq=np.r_[q_joint,1.],bounds=(0,None),method='highs')
    tests['inside_joint_training_feature_hull']=bool(joint.success)
    extra=not all(tests.values())
    if extra and not allow_extrapolation:
        return {'state':'OUT_OF_TRAINING_SUPPORT','query':row,'support_checks':tests,
                'message':'Use --allow-extrapolation only for explicitly exploratory predictions; no optimization-grade validity.'}
    ray=ray_predict(row,rays);q=dict(row,ray_A_total=ray['A_total'])
    m=joblib.load(HERE/f'models/{pool}__pooled__{method}.joblib')
    p,s=predict(m,[q])
    value=float(p[0])
    return {'state':'EXPLORATORY_EXTRAPOLATION' if extra else 'EXPLORATORY_INTERPOLATION',
       'query':row,'support_checks':tests,'method':method,'pool':pool,'predicted_A_loss':value,
       'GP_std_uncalibrated':float(s[0]) if s is not None else None,
       'ray_A_total':ray['A_total'],'ray_A_film':ray['A_film'],'ray_A_Ti_total':ray['A_Ti_total'],
       'ray_sampling_standard_error':ray['se_A_total'],'physical_range_0_1':0<=value<=1,
       'elapsed_s':time.monotonic()-start,
       'meaning':'Single wavelength/direction/polarization absorptivity; directional emissivity under reciprocity and thermal-equilibrium assumptions.',
       'limitations':['Not a band-integrated or cone-integrated result.','Eight training geometries; changed bottom width unsupported.',
          'No clipping of unphysical predictions. GP uncertainty and ray sampling error do not measure model discrepancy.',
          'Normal-angle wave effects and unreviewed training labels remain limitations.']}

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--height',type=float,required=True);p.add_argument('--film',type=float,required=True)
    p.add_argument('--bottom',type=float,default=15);p.add_argument('--wavelength',type=float,default=10.5)
    p.add_argument('--theta',type=float,default=0);p.add_argument('--phi',type=float,default=0)
    p.add_argument('--pol',choices=['s','p'],default='s');p.add_argument('--allow-extrapolation',action='store_true')
    p.add_argument('--method',default='gp');p.add_argument('--pool',default='primary_exploratory')
    a=p.parse_args();q=dict(height_um=a.height,film_um=a.film,bottom_width_um=a.bottom,
       wavelength_um=a.wavelength,theta_deg=a.theta,phi_deg=a.phi,polarization=a.pol)
    print(json.dumps(evaluate(q,a.allow_extrapolation,a.method,a.pool),ensure_ascii=False,indent=2,allow_nan=False))

if __name__=='__main__': main()
