# Tail diagnostic offline report recovery

Diagnostic only. Original measurements, original status and original reports are unchanged.
Recovered FFT uses verified regular samples; no interpolation or resampling.

Original run state: DIAGNOSTIC_COMPLETE.
Original H qualification: TIME_LIMIT_UNQUALIFIED.
Control reproduction gate passed: True.

- pml6_control: DIAGNOSTIC_COMPLETE; qualification=false.
  Authoritative supervisor status=DIAGNOSTIC_COMPLETE; worker status=DIAGNOSTIC_COMPLETE.
  Regular samples: 4800; original samples: 4800.
- pml12_test: DIAGNOSTIC_COMPLETE; qualification=false.
  Authoritative supervisor status=DIAGNOSTIC_COMPLETE; worker status=DIAGNOSTIC_COMPLETE.
  Regular samples: 4800; original samples: 4800.

Matched H through t=1200.0; required t>=1000; maximum normalized error=2.6637531547534046e-15.

Exact matched PML6/PML12 interval: [1000.0, 1200.0]; 801 samples.
Mean back intensity / H reference: control=0.005603725878299913; thick PML=0.0034733174741095655; thick/control=0.6198228731279979.
Maximum back intensity / H reference: control=0.006093013986960662; thick PML=0.0037178012376565804; thick/control=0.6101744137815621.
Matched-window signed complex spectra are recorded in diagnostic_summary.json.
PML-thickness sensitivity diagnostic only. Differences do not by themselves establish causality, convergence, or optical qualification.

These diagnostics do not establish a PML-boundary cause or a qualified optical size trend.
Details and evidence hashes: diagnostic_summary.json.
