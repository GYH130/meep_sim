# 离散平面 TM 反射诊断的可复核推导

此文件只解释保存的平面 r16 结果，不是新 FDTD、斜槽验证或误差校正方法。
长度单位 μm，角频率 ω=2π/λ；向金属内的深度为 z，Ex 在 z=jh，Hz 在 z=(j+1/2)h。
对 exp(−iωt) 的 Yee 更新，定义 Ω=2 sin(ωΔt/2)/Δt，K=2 sin(kx h/2)/h。

对每个 Lorentz/Drude 极点，中心差分极化更新给出

εd = ε∞ + Σ σω₀²/[ωden² − Ω² − iγ sin(ωΔt)/Δt]。

其中所有 ω₀、γ 都先由输入循环频率乘 2π；Lorentz 的 ωden=ω₀，Drude 的 ωden=0。
此式使用冻结材料参数，未从保存的 D/E 或 R 拟合。

令 Bm=Ω−K²/(Ωεd)，Ba=Ω−K²/Ω。在均匀金属/空气中的每格传播因子分别为

p=exp[2i asin(h√(Ω²εd−K²)/2)]，a=exp[2i asin(h√(Ω²−K²)/2)]。

被动半空间选择 |p|≤1 的支。消去 Hz 后，对 Ex 的二阶递推为
E(j+1)−2E(j)+E(j−1)+h²(Ω²εd−K²)E(j)=0，因此 p+1/p−2=−h²(Ω²εd−K²)。

当前平面对齐在 Ex 节点；只平滑瞬时 ε∞、不平滑色散项时，界面节点的近似构成值为
εI=εd−(ε∞−1)/2。保存的界面 D/E 可作单独辅助回放，但不用于默认独立预测。

界面左侧为空气、右侧为金属；节点处的离散更新为

(E1−E0)/Bm − (E0−E−1)/Ba + Ωh²εI E0 = 0。

代入 E0=1+r，E1=p(1+r)，E−1=a⁻¹+ra，令 L=(p−1)/Bm+Ωh²εI，得到

r=[1−a⁻¹−BaL]/[BaL−(1−a)]，R=|r|²。

自检：空气无反差时 r=0；h→0 时趋于连续 TM Fresnel；Δt→0 时 εd 趋于原材料模型。
实际膜厚 60 μm 远大于约 0.081 μm 的振幅衰减长度，本诊断采用半空间近似；正式验收仍对完整 FDTD 与连续解析膜解比较。

限制：界面瞬时平滑采用理想半占据；实际保存的界面值存在很小差别，因此不声称逐位复现。
预测恰好符合某个平面门槛不能证明斜槽、偏振、镜像或网格收敛通过，不能把预测误差从结果中扣除来放行。

依据：[Meep 色散界面处理](https://meep.readthedocs.io/en/latest/Subpixel_Smoothing/#what-about-dispersive-materials)，[中心差分极化更新](https://github.com/NanoComp/meep/blob/master/src/susceptibility.cpp)。
