"""Local exploratory forward surrogate; never connects to a simulation server."""
import os
os.environ.setdefault('OMP_NUM_THREADS', '1')
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
from pathlib import Path
import json, warnings
import numpy as np
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.dummy import DummyRegressor
from sklearn.linear_model import Ridge
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel, Matern, WhiteKernel

FEATURE_NAMES = ['height_um', 'bottom_width_um', 'film_um', 'wavelength_um',
                 'cos_theta', 'sin_theta_cos_6phi', 'is_p']
METHODS = ['mean', 'ridge', 'gp', 'ray', 'ray_mean_residual', 'ray_ridge_residual', 'ray_gp_residual']

def features(rows, mode='pooled'):
    if mode == 'condition':
        return np.array([[r['height_um'], r['film_um']] for r in rows], float)
    return np.array([[r['height_um'], r['bottom_width_um'], r['film_um'], r['wavelength_um'],
                      np.cos(np.deg2rad(r['theta_deg'])),
                      np.sin(np.deg2rad(r['theta_deg']))*np.cos(6*np.deg2rad(r['phi_deg'])),
                      float(r['polarization']=='p')] for r in rows], float)

def condition(r):
    return f"l{r['wavelength_um']:g}_th{r['theta_deg']:g}_phi{r['phi_deg']:g}_{r['polarization']}"

def fit(rows, method, mode='pooled'):
    out = dict(method=method, mode=mode, train_designs=sorted({r['design_id'] for r in rows}),
               n_rows=len(rows), feature_names=FEATURE_NAMES if mode=='pooled' else FEATURE_NAMES[::2][:2],
               status='EXPLORATORY_NOT_FORMALLY_QUALIFIED', warnings=[])
    if method == 'ray':
        out['model'] = None
        return out
    base = method.split('_')[1] if method.startswith('ray_') else method
    if base == 'mean':
        estimator = DummyRegressor(strategy='mean')
    elif base == 'ridge':
        estimator = Ridge(alpha=.1)
    else:
        kernel=ConstantKernel(1., (1e-3,1e3))*Matern(1.,(.05,20.),nu=2.5)+WhiteKernel(.01,(1e-5,1.))
        estimator=GaussianProcessRegressor(kernel=kernel,normalize_y=True,alpha=1e-10,
                        random_state=20261102,n_restarts_optimizer=2 if mode=='condition' else 0)
    y=np.array([r['A_loss_total'] for r in rows])
    if method.startswith('ray_'):
        y-=np.array([r['ray_A_total'] for r in rows])
    model=make_pipeline(StandardScaler(),estimator)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always');model.fit(features(rows,mode),y)
    out.update(model=model,warnings=sorted({str(w.message) for w in caught}))
    return out

def predict(model, rows):
    if model['method']=='ray':
        return np.array([r['ray_A_total'] for r in rows]),None
    x=features(rows,model['mode'])
    if 'gp' in model['method']:
        p,s=model['model'].predict(x,return_std=True)
    else:
        p=model['model'].predict(x);s=None
    if model['method'].startswith('ray_'):
        p=p+np.array([r['ray_A_total'] for r in rows])
    return p,s

def metrics(records):
    e=np.array([r['predicted']-r['observed'] for r in records])
    y=np.array([r['observed'] for r in records])
    p=np.array([r['predicted'] for r in records])
    return dict(n=len(e),n_designs=len({r['heldout_design'] for r in records}),
                MAE_pp=float(abs(e).mean()*100),RMSE_pp=float(np.sqrt((e*e).mean())*100),
                max_abs_error_pp=float(abs(e).max()*100),bias_pp=float(e.mean()*100),
                R2=float(1-np.sum(e*e)/np.sum((y-y.mean())**2)) if len(y)>1 and np.std(y)>0 else None)
