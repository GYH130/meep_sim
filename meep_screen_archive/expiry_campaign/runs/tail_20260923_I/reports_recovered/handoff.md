# Tail diagnostic offline report recovery

Diagnostic only. Original measurements, original status and original reports are unchanged.
Recovered FFT uses verified regular samples; no interpolation or resampling.

Original run state: FAILED.
Original H qualification: TIME_LIMIT_UNQUALIFIED.
Control reproduction gate passed: False.

- pml6_control: TIME_LIMIT_UNQUALIFIED; qualification=false.
  Authoritative supervisor status=TIME_LIMIT_UNQUALIFIED; worker status=TIME_LIMIT_UNQUALIFIED.
  Regular samples: 3658; original samples: 3659.
  FFT-only exclusion: terminal t=914.625, interval=0.125 instead of 0.25; original NPZ retained.
- pml12_test: NOT_RUN_OR_NO_RESULT; qualification=false.

Matched H through t=914.0; required t>=1000; maximum normalized error=2.6637531547534046e-15.

NOT_EVALUATED: requires completed control and PML12 arms plus original control reproduction gate

These diagnostics do not establish a PML-boundary cause or a qualified optical size trend.
Details and evidence hashes: diagnostic_summary.json.
