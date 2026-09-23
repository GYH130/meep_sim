"""Presentation-only postprocessor: show post-source tails on readable axes."""
from pathlib import Path
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    run=Path(__file__).resolve().parent/'runs/tail_20260923_I'
    old=json.loads((run/'inputs/H_result.json').read_text())
    trace=json.loads((run/'inputs/H_stop.json').read_text())['trace']
    source_end=old['source_end_time']
    tail=[r for r in trace if r['time']>=source_end]
    floor=old['stop']['normalization']['denominator_floor']
    fig,axs=plt.subplots(1,2,figsize=(12,4.3),constrained_layout=True)
    for n in ('front_air','back_air','cavity_0.5'):
        axs[0].semilogy([r['time'] for r in tail],[r['probe_ratio'][n] for r in tail],label=n)
    axs[0].set(title='H: post-source tail (unqualified)',xlabel='Simulation time (um/c)',ylabel='Intensity / stop denominator',ylim=(1e-7,.1))
    axs[0].legend(fontsize=8)
    plotted=False
    for arm in ('pml6_control','pml12_test'):
        path=run/arm/'power_trace.json'
        if path.exists():
            rows=[r for r in json.loads(path.read_text()) if r['time']>=source_end]
            if rows:
                plotted=True
                axs[1].semilogy([r['time'] for r in rows],[r['probe_power']['back_air']/floor for r in rows],label=arm)
    axs[1].set(title='I: PML thickness comparison (diagnostic)',xlabel='Simulation time (um/c)',ylabel='Back |E|² / H excitation peak',ylim=(1e-7,.1),xlim=(source_end,1200))
    if plotted: axs[1].legend(fontsize=8)
    else: axs[1].text(.5,.5,'Post-source diagnostic pending',ha='center',transform=axs[1].transAxes)
    for ax in axs:
        ax.axhline(1e-4,color='black',ls='--',linewidth=1)
        ax.grid(True,alpha=.2)
    fig.savefig(run/'reports/tail_comparison.png',dpi=160)
    plt.close(fig)


if __name__=='__main__': main()
