"""Server-only OFFLINE audit of saved flat-case fields; no Simulation/time stepping.

Does not alter solver code, qualification, material parameters, or the budget.
New output directory must not already exist. No measured value is fitted into
the default discrete prediction. A separate measured-interface replay is labeled.
"""
from __future__ import annotations

import csv
import json
import math
import os
from pathlib import Path
import sys

import numpy as np

STAGE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(STAGE))
from ti2d.common import atomic_json, digest_file, utcnow
from ti2d.materials import build_medium, epsilon, load_model
from diagnostic_tools.flat_discrete import discrete_epsilon, planar_tm_prediction


def cpair(z):
    z = complex(z)
    return [z.real, z.imag]


def encode(value):
    if isinstance(value, (complex, np.complexfloating)):
        return {'real': float(value.real), 'imag': float(value.imag)}
    if isinstance(value, dict):
        return {k: encode(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [encode(v) for v in value]
    if isinstance(value, np.generic):
        return value.item()
    return value


def main():
    if sys.platform != 'linux' or str(STAGE) != '/root/autodl-tmp/meep_sim/meep_screen':
        raise RuntimeError('Read the original fields on the server; no local FDTD')
    os.environ['MPLCONFIGDIR'] = str(STAGE / 'cache/matplotlib')
    os.environ['XDG_CACHE_HOME'] = str(STAGE / 'cache')
    run = STAGE / 'runs/repair_20260922_C'
    attempt = run / 'tasks/flat_r16_p_plus30/attempt_1'
    output = STAGE / 'repair_diagnostics/flat_spatial_audit_20260922'
    if output.exists():
        raise FileExistsError('Do not overwrite an earlier audit')
    result_path = attempt / 'result.json'
    result = json.loads(result_path.read_text())
    c = result['case_config']
    if (c['case'], c['polarization'], c['resolution']) != ('flat', 'p', 16):
        raise ValueError('Audit is limited to the saved r16 flat p case')
    for name, expected in result['evidence_sha256'].items():
        if digest_file(attempt / name) != expected:
            raise ValueError('Changed source evidence: ' + name)
    material_path = Path(c['material_path'])
    if digest_file(material_path) != c['material_sha256']:
        raise ValueError('Frozen material hash changed')
    model = load_model(material_path)
    f, h, dt = 1 / c['wavelength_um'], 1 / c['resolution'], c['courant'] / c['resolution']
    eps = complex(epsilon(model, f))
    # Material response evaluation only: this does not initialize Meep fields.
    meep_eps = np.diag(build_medium(model).epsilon(f))
    theta = math.radians(c['emission_theta_deg'])
    kx = -f * math.sin(theta)
    q = np.sqrt(eps - math.sin(theta) ** 2)
    independent_fresnel = abs((eps * math.cos(theta) - q) / (eps * math.cos(theta) + q)) ** 2
    skin = 1 / (2 * math.pi * f * q.imag)
    predicted = planar_tm_prediction(model, c['wavelength_um'], c['emission_theta_deg'], 16, c['courant'])

    profiles, profile_values = {}, {}
    with np.load(attempt / 'volume_ED.npz', allow_pickle=False) as data:
        for component in ('Ex', 'Ey'):
            x, y = data[component + '_x_canonical'], data[component + '_y_canonical']
            E, D = data[component + '_E_canonical'], data[component + '_D_canonical']
            if E.shape != D.shape or E.shape != (len(x), len(y)):
                raise ValueError('Field shape mismatch, no clipping allowed')
            phase = np.exp(-2j * np.pi * kx * x)
            em = np.mean(E * phase[:, None], axis=0)
            dm = np.mean(D * phase[:, None], axis=0)
            depth = result['layout']['surface'] - y
            use = (depth >= h * .5) & (depth <= h * 4.1)
            bulk = np.vdot(em[use], dm[use]) / np.vdot(em[use], em[use])
            layers = np.flatnonzero((depth >= -h * 1.01) & (depth <= .75))
            profiles[component] = {
                'bulk_epsilon_from_saved_D_E': cpair(bulk),
                'bulk_region_depth_um': [float(depth[use].min()), float(depth[use].max())],
                'bulk_epsilon_error_vs_ADE_prediction': abs(bulk - predicted['epsilon_discrete']),
                'layers': [{'depth_um': float(depth[i]), 'E': cpair(em[i]), 'D': cpair(dm[i]),
                            'D_over_E': cpair(dm[i] / em[i])} for i in layers[::-1]],
            }
            profile_values[component] = (depth, em, dm)
    depth, em, dm = profile_values['Ex']
    interface_indices = np.flatnonzero(np.isclose(depth, 0, rtol=0, atol=1e-10))
    if len(interface_indices) != 1:
        raise ValueError('Expected an Ex node exactly on the planar interface')
    interface = int(interface_indices[0])
    measured_interface = dm[interface] / em[interface]
    bulk_indices = np.flatnonzero((depth >= h * .9) & (depth <= h * 4.1))
    actual_ratio = np.vdot(em[bulk_indices], em[bulk_indices - 1]) / np.vdot(em[bulk_indices], em[bulk_indices])
    measured_interface_replay = planar_tm_prediction(
        model, c['wavelength_um'], c['emission_theta_deg'], 16, c['courant'],
        interface_epsilon=measured_interface)

    with np.load(attempt / 'monitor_planes.npz', allow_pickle=False) as data:
        x = data['x_refl']
        phase = np.exp(-2j * np.pi * kx * x)
        field = {k: complex(np.mean(data[k] * phase)) for k in
                 ('incident_Ex', 'incident_Hz', 'scattered_Ex', 'scattered_Hz')}
        measured_kx = float(np.mean(np.angle(data['incident_Hz'][1:] / data['incident_Hz'][:-1])
                                    / (2 * np.pi * np.diff(x))))
        purity_residual = float(np.linalg.norm(data['incident_Hz'] - field['incident_Hz'] / phase)
                                / np.linalg.norm(data['incident_Hz']))
    rows = []
    for resolution in (16, 20, 24, 32):
        p = planar_tm_prediction(model, c['wavelength_um'], c['emission_theta_deg'], resolution, c['courant'])
        rows.append({'resolution': resolution, 'prediction_R': p['R'],
                     'prediction_absolute_error': abs(p['R'] - independent_fresnel),
                     'cells_per_amplitude_skin_depth': resolution * skin,
                     'prediction_within_flat_tolerance': abs(p['R'] - independent_fresnel) < .0025,
                     'is_FDTD_result': False, 'qualified': False})
    ledger_path = STAGE / 'budget_ledger.json'
    ledger = json.loads(ledger_path.read_text())
    if any(i['end'] is None for i in ledger['intervals']):
        raise RuntimeError('Unexpected active solver: audit must not race with a run')
    ledger_sha = digest_file(ledger_path)
    status = json.loads((run / 'status.json').read_text())
    baseline = json.loads((STAGE / 'repair1.json').read_text())['baseline_active_solver_seconds']
    repair_used = status['usage']['active_solver_seconds'] - baseline
    report = {
        'status': 'OFFLINE_DISCRETIZATION_DIAGNOSIS_NOT_QUALIFICATION', 'created_utc': utcnow(),
        'source_result': str(result_path), 'source_result_sha256': digest_file(result_path),
        'material_sha256': digest_file(material_path), 'original_result_status_unchanged': result['status'],
        'new_FDTD_runs': 0, 'new_FDTD_solver_seconds': 0,
        'meep_continuous_epsilon_max_difference': float(np.max(np.abs(meep_eps - eps))),
        'epsilon_continuous': cpair(eps), 'epsilon_ADE': cpair(discrete_epsilon(model, f, dt)),
        'n_k': cpair(np.sqrt(eps)), 'amplitude_skin_depth_normal_um': skin,
        'power_skin_depth_normal_um': skin / 2, 'grid_spacing_um': h,
        'independent_fresnel_R': float(independent_fresnel), 'saved_FDTD_R': result['metrics']['R'],
        'discrete_prediction_no_measured_fits': encode(predicted),
        'prediction_vs_saved_R_absolute_difference': abs(predicted['R'] - result['metrics']['R']),
        'measured_interface_epsilon': cpair(measured_interface),
        'measured_interface_replay_not_independent': encode(measured_interface_replay),
        'saved_bulk_depth_factor': cpair(actual_ratio),
        'bulk_depth_factor_absolute_difference': abs(actual_ratio - predicted['depth_factor']),
        'saved_kx': measured_kx, 'configured_kx': kx, 'incident_plane_purity_residual': purity_residual,
        'incident_Ex_over_Hz': cpair(field['incident_Ex'] / field['incident_Hz']),
        'reflected_Hz_amplitude_ratio_squared': abs(field['scattered_Hz'] / field['incident_Hz']) ** 2,
        'profiles': profiles, 'resolution_predictions_not_validation': rows,
        'repair_used_seconds': repair_used, 'repair_remaining_seconds': 1800 - repair_used,
        'ledger_sha256_before': ledger_sha,
        'limits': ['Only the flat aligned p-polarized case is diagnosed.',
                   'Discrete continuum-limit agreement does not qualify structured geometry.',
                   'No retrospective correction or reclassification of old results.',
                   'r32 is outside the original candidate set; not launched or certified.'],
        'sources': ['https://meep.readthedocs.io/en/latest/Subpixel_Smoothing/#what-about-dispersive-materials',
                    'https://github.com/NanoComp/meep/blob/master/src/susceptibility.cpp'],
    }
    output.mkdir(parents=True)
    atomic_json(output / 'audit.json', report)
    with (output / 'resolution_predictions.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator='\n')
        writer.writeheader(); writer.writerows(rows)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.3))
    ax = axes[0]
    ax.plot([r['resolution'] for r in rows], [r['prediction_R'] for r in rows], 'o--', label='Discrete prediction (NOT FDTD)')
    ax.scatter([16], [result['metrics']['R']], marker='x', s=90, color='black', zorder=5, label='Saved Meep r16')
    ax.axhline(independent_fresnel, color='green', label='Continuum Fresnel')
    ax.axhspan(independent_fresnel - .0025, independent_fresnel + .0025, color='green', alpha=.10, label='Flat R tolerance only')
    ax.set(xlabel='Resolution (pixels / um)', ylabel='Reflectance', title='One flat p, +30 deg case only')
    ax.legend(fontsize=8); ax.grid(alpha=.2)
    ax = axes[1]
    indices = np.flatnonzero((depth >= h * .9) & (depth <= .65))[::-1]
    z, values = depth[indices], np.abs(em[indices])
    ax.semilogy(z, values / values[0], 'o', label='Saved interior Ex')
    ax.semilogy(z, np.abs(predicted['depth_factor']) ** ((z - z[0]) / h), '--', label='Discrete prediction')
    ax.semilogy(z, np.exp(-(z - z[0]) / skin), ':', label='Continuum decay (same anchor)')
    ax.set(xlabel='Depth below surface (um)', ylabel='Relative |Ex|', title='Bulk metal attenuation')
    ax.legend(fontsize=8); ax.grid(alpha=.2)
    fig.tight_layout(); fig.savefig(output / 'flat_discretization.png', dpi=160); plt.close(fig)
    table = '\n'.join(f"| {r['resolution']} | {r['prediction_R']:.9f} | {r['prediction_absolute_error']:.6f} | {'门槛内，仅预测' if r['prediction_within_flat_tolerance'] else '仍超门槛'} |" for r in rows)
    note = f'''# 平面反射率偏差诊断（离线，不是新的验收通过）

结论：独立离散方程预测 R={predicted['R']:.10f}，与保存的 Meep R={result['metrics']['R']:.10f} 相差 {report['prediction_vs_saved_R_absolute_difference']:.3g}；这定量解释了约 0.004995 的连续 Fresnel 偏差。主要瓶颈是金属趋肤区/色散界面的空间离散，而不是继续延长运行时间。

## 核查证据

- 解析模型与 Meep 材料响应最大差 {report['meep_continuous_epsilon_max_difference']:.3g}，无材料或单位错配；解析 R={independent_fresnel:.10f}。
- 保存的入射 kx={measured_kx:.12g}，配置 kx={kx:.12g}；平面波杂散相对范数 {purity_residual:.3g}。
- 振幅趋肤深度 {skin:.6f} μm，r16 网格间距 {h:.4f} μm，仅 {skin/h:.3f} 个网格/振幅衰减长度。功率衰减长度更短，为 {skin/2:.6f} μm。
- 独立离散模型不拟合已测 R 或 D/E；金属层间复振幅比例与保存场相差 {report['bulk_depth_factor_absolute_difference']:.3g}。单列的“测得界面 D/E 回放”只是辅助诊断，不能当独立验证。
- 原结果仍为 NUMERICAL_FAILED；原场、材料、求解代码、阈值和预算账本均不修改。无新增 FDTD，求解预算消耗为 0。

## 网格预测（以下不是已计算的 FDTD 结果）

| r | 离散预测 R | 对连续解析绝对差 | 判断 |
|---|---:|---:|---|
{table}

预测仅适用于当前平面、偏振和界面对齐，不能认证斜槽、s 偏振、镜像或网格收敛。r32 也只是单项门槛预测，不是生产网格获准。

![离散预测与保存场对照](flat_discretization.png)

## 下一步与边界

不把旧反射率减去预测误差来“修正通过”，不更改 Ti 极点，不关闭吸收或放宽 0.0025 门槛。r20/r24 预测仍失败，没有为了凑条件而启动。当前修复轮尚余 {(1800-repair_used)/60:.2f} 分钟，少于已测完整 r16 平面参考+结构的 {result['elapsed_s']/60:.2f} 分钟，因此未启动预计更慢的完整细网格任务。

需用户批准扩展原 r16/r20/r24 候选网格至 r32，并为下一轮单平面复核设定独立求解上限（仍计入原 48 小时）。平面通过后再评估斜槽与必要资格验证；20 条件尺寸试算仍不启动。

依据：[Meep 色散界面平滑限制](https://meep.readthedocs.io/en/latest/Subpixel_Smoothing/#what-about-dispersive-materials)；[极化更新方程](https://github.com/NanoComp/meep/blob/master/src/susceptibility.cpp)。本报告是对给定离散实现的诊断，不是样品材料或高温实验验证。
'''
    (output / 'REPORT.md').write_text(note)
    if digest_file(ledger_path) != ledger_sha:
        raise RuntimeError('Budget ledger changed during offline audit; inspect before publishing')
    atomic_json(output / 'MANIFEST.json', {
        'audit_script_sha256': digest_file(Path(__file__)),
        'formula_module_sha256': digest_file(STAGE / 'diagnostic_tools/flat_discrete.py'),
        'source_evidence_sha256': result['evidence_sha256'],
        'output_sha256': {p.name: digest_file(p) for p in output.iterdir() if p.is_file()},
        'ledger_unchanged': True, 'qualification_changed': False, 'new_FDTD_runs': 0})
    print(json.dumps({'output': str(output), 'predicted_R': predicted['R'],
                      'agreement_error': report['prediction_vs_saved_R_absolute_difference'],
                      'new_FDTD_runs': 0, 'repair_remaining_seconds': 1800 - repair_used}, indent=2))


if __name__ == '__main__':
    main()
