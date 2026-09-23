# 2026-09-23 尾场诊断 I

目的：定位 H 轮斜槽基准的背面慢衰减，不把超时数据升级为合格结果。

- 位置：服务器 `/root/autodl-tmp/meep_sim/meep_screen/expiry_campaign/runs/tail_20260923_I`。
- 原状态：H 平面 r32 p 条件通过；斜槽 r32 p 六小时超时，背面末端强度比约 0.00261，阈值 0.0001。
- 第 1 组：PML=6 μm，复现原斜槽轨迹。被动 DFT/体场累加器省略；几何、材料、光源和离散参数一致。必须与 H 原始多探针轨迹匹配，才允许执行第 2 组。
- 第 2 组：PML=12 μm。非 PML 区域、光源、槽形及探针物理坐标保持不变。
- 两组：Ordal Route A、r32、Courant=0.25、eps_averaging=false、p 偏振、观察角 +30°；每组 4 MPI rank、单线程、顺序运行；均从零场开始。
- 每组目标 t=1200（源 t=525 结束），墙钟上限 2.5 小时；合计最多 5 小时，并计入原 48 小时累计求解账本。
- 2026-09-23 用户更正租期至明天，撤销今天 16:52:44 的调度截止。见 `RENTAL_DEADLINE_OVERRIDE.json`。准确到期时刻尚未提供。
- 当前计算已获得完整 2.5 小时单次额度，正在运行的进程按已分配时长执行，不需要重启。冻结执行清单保留旧截止用于追溯，后续调度以更正文件为准。本轮 5 小时和总计 48 小时求解预算不变。
- 输出：9 点 Ex/Ey/Hz 复波形（dt=0.25）、3 条横向复场线（dt=2）、非 PML 峰值位置（dt=20）、少量场强图（dt=100）。保存原始数组、坐标、时间及 SHA256。
- 结果只能说明尾场频率、空间分布及对吸收边界的敏感性。本轮不生成正式 R/T/A 或尺寸敏感性。
- 测试：四项离线检查覆盖频谱、历史轨迹复现判据、预算截止及物理探针坐标一致性。

完成文件：`reports/diagnostic_summary.json`、`reports/handoff.md`、`reports/tail_comparison.png`、`reports/tail_spectra.png`。
短状态：`status.json`；当前算例进度：`pml6_control/heartbeat.json` 或 `pml12_test/heartbeat.json`。

后台队列完成后生成 `expiry_campaign/exports/I_raw_evidence.tar.gz`。本机单次归档程序将下载并校验该文件，再向已授权的 `codex/archive-20260922-ti2d` 归档分支追加脱敏代码、配置、小结果和报告；不会升级科学资格。
本机需保持开机、联网，直至 `meep_screen/backups/20260922_expiry/I_tail_diagnostic/collector_status.json` 为 COMPLETE。计算本身独立运行于服务器。
后台求解、分析、归档程序不调用 AI API。未取得实际账单，不估造 token 数或人民币成本。
