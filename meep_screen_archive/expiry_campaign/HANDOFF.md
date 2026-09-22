# Expiry validation H — 2026-09-22

## Decision and boundaries

User authorizes more bounded server work and a NEW sanitized GitHub archive branch. No merge, force push, modifications of old results or main. User-reported rental remaining 24h from 2026-09-22 12:52:44 UTC; not independently confirmed with provider.

Core fork lives ONLY here. Original `../ti2d` remains unchanged (24d09c50dcf6c0cb9271947dbcb04f03b62d9d0500bd7bace34e6909d54834c1). Material remains Ordal Route A (4791b3868810fb4a3fa93c8469b488650270c0f47a18804c9b26bd7354153b5e). The sole physical/numerical option change is explicit `eps_averaging=False`. All incident source, geometry, E/D integration, order projection and stop algorithms retained. Cache identity includes the new implementation and smoothing option; no old-reference reuse.

Old D flat r32 p+30 is individually qualified. E slanted r32 overflowed. F/G short tests suggest smoothing-off removes observed EARLY growth, but G stopped at t98.83 before source end525: full stability NOT established. Old results and statuses stay unchanged.

## Locked work

Ten sequential validation conditions: flat32p, slot32p, strict32p (1e-6), mirror32p (tilt−30/observation−30), flat32s, slot32s, flat24p, slot24p, flat20p, slot20p. All at 10.5µm, C0.25, 4 ranks ×1 thread. Each depends on the preceding qualified result. Fixed individual scientific gates; strict/mirror/mesh comparisons checked as available. Any failed required case/comparison halts expansion, no automatic physical-failure retry.

This is ONLY a qualification subset; production remains locked even if all ten pass. Straight-slot checks, full s symmetry, L48 boundary and production-grid evidence remain missing. No 20-condition size study or 90-condition expansion this turn. r32 is not declared production-qualified.

## Runtime

Original shared 48h union ledger remains authoritative (before H:32365.66657280922s); H additionally limited to20h and absolute 2026-09-23 08:52:44 UTC /16:52:44 China. Admission includes remaining deadline and30s cleanup reserve. Session guard signals only own subprocess identities and verifies orphan MPI ranks are gone before closing ledger. Memory≤70% current container available; never host free-memory values.

No field checkpoint resume. Only hash-complete qualified cache may be reused. Interrupted/unqualified previous attempts require human review rather than silent restart. No AI/API calls from the queue; no inferred RMB cost. Four-hour rental reserve covers packaging and off-server backup, not more simulation.

## Commands (server only)

Existing environment `/root/miniconda3/envs/meep_sim`, no install/upgrade.

```sh
cd /root/autodl-tmp/meep_sim/meep_screen/expiry_campaign
/root/miniconda3/envs/meep_sim/bin/python -B create_campaign.py
nohup /root/miniconda3/envs/meep_sim/bin/python -B run_campaign.py > supervisor.log 2>&1 < /dev/null &
```

Create ONCE only. Do not run another supervisor. Short status:

```sh
cd /root/autodl-tmp/meep_sim/meep_screen/expiry_campaign
/root/miniconda3/envs/meep_sim/bin/python -c 'import json; p=json.load(open("runs/expiry_20260922_H/status.json")); print({k:p.get(k) for k in ("state","active_task","remaining_solver_seconds","tasks")})'
```

Queue state: `runs/expiry_20260922_H/status.json`; small reports in that run's reports directory. Per-condition results and raw arrays under `tasks/<id>/attempt_1`. `exports/compact_latest.tar.gz` rebuilt after each condition and queue exit. `exports/H_raw_evidence.tar.gz` +SHA256 generated after supervisor exit. Check `exports/final_status.json` for packaging failures. Packaging/backup is separate from the solver budget.

## Backup and GitHub

Before H, server archive `../exports/preexpiry_evidence_20260922.tar.gz` copied to local backup directory and SHA verified before claiming success. Compact export scans token/private-key patterns and verifies qualified raw-evidence hashes. It excludes full NPZ/HDF5 arrays, caches, credentials and full logs. Preserve the separate RAW archive off-server; GitHub small report alone is NOT a full field backup.

Original dirty checkout is never pulled, merged, switched or committed by this campaign. Use isolated local archive checkout, new branch `codex/archive-20260922-ti2d`; push only with verified authentication and no force. Follow-up retrieves new H results before rental expiry, without waiting/polling through a long model conversation.
