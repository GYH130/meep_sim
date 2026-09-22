# E: one baseline slanted-slot validation, not a production scan

## Frozen scope

Run: `baseline_r32_20260922_E`; task: `baseline_r32_p_plus30`.
Server root: `/root/autodl-tmp/meep_sim/meep_screen`.
Only the flat/structure switch changes the physical problem relative to the
qualified `flat_r32_20260922_D` case. The task ID and approved wall-time cap also
change. No edits to the `ti2d` solver, frozen Ordal Route A material, acceptance
thresholds, MPI layout, source, or stopping policy are authorized by this run.

- P=50 um, surface fill=0.692820323, axial length=40 um, tilt=30 degrees;
  Ti thickness=60 um; no oxide layer.
- Resolution=32 pixels/um, Courant=0.25, wavelength=10.5 um, p polarization,
  emission observation angle=+30 degrees. Incident k is stored separately.
- One job, four MPI ranks, one thread per rank. No automatic retry or follow-on.
- An exactly matched, checksum-qualified incident reference is available from
  run D. Reuse is not reuse of the flat target result; the slanted structure
  must be solved afresh. Any fresh reference work also counts against the cap.

## Admission evidence and budget

The r32 flat case converged and passed all single-case gates: Fresnel R error
0.0020983751 < 0.0025; independent absorption difference 9.6224e-8; reflection
order closure error 2.3558e-5. These do not certify slot or production accuracy.

Shared historical solver time at admission: 30121.365798711777 seconds.
This single task has at most 21600 seconds (six hours), including reference,
target, failures, and shutdown; effective shared ceiling 51721.36579871178
seconds, also inside the original 172800-second stage ceiling. History is not
reset. The prior run D remains a separate three-hour-limited phase.
Use container CPU/memory limits, not host `free` output. The process-tree memory
limit is 70% of available container memory at launch.

## What to inspect after completion

`runs/baseline_r32_20260922_E/tasks/baseline_r32_p_plus30/attempt_1/result.json`
must show a converged stop, independent absorption/flux consistency, upward
reflection and downward transmission order closure, finite passive metrics,
and qualified background reference. Preserve raw signed values and all failed
evidence. TIME_LIMIT_UNQUALIFIED is never a pass.

Even if this task is QUALIFIED, the campaign gate remains incomplete until the
locked mesh, mirror, polarization, and strict-stop evidence is complete. No
20-condition pilot, feature-size sensitivity, production grid qualification,
or validated uncertainty bound is claimed by this run alone.

## Operation and outputs

Preflight (no FDTD):
`python -B -m ti2d.cli preflight --manifest runs/baseline_r32_20260922_E/manifest.json`

Launch once in a detached server process:
`python -B -m ti2d.queue --manifest runs/baseline_r32_20260922_E/manifest.json --phase validation`

Short status:
`python -B -m ti2d.cli status --manifest runs/baseline_r32_20260922_E/manifest.json`

Automatic report: `runs/baseline_r32_20260922_E/reports/handoff.md`.
The generic report describes campaign-level gates; inspect the individual task
status separately. `create_baseline_r32_run.py` refuses overwrite, missing or
altered history, changed source physics/material/implementation, and incomplete
source evidence. An existing unqualified attempt requires review, not resume.

AI ends the turn after one startup verification. No AI polling, API requests
from the queue, automatic desktop notification, Git push, or old-file deletion.
The queue persists across SSH disconnects. Actual AI billing is not available
from the solver and must not be inferred in RMB.
