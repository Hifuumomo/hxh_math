# 仅水平黏性可压缩流体研究审计

评估日期：2026 年 10 月 30 日。对象是当前 `research/` 中的研究源稿及验证脚本。

**未发现足以推翻指定测试初值的二阶声学生成、非抵消公式或对数时间交接论证的明确错误。但当前项目尚未证明该真正二维初值会形成激波，也未证明所猜想的寿命渐近常数。** 局部能量部分是一份可以补全的证明草图；长期调制和奇性形成则仍有实质性证明工作。稿件已经明确区分这些层次，不能把明确标注的猜想误报成已被反例推翻的定理。

本次保留原研究文件，新增独立检查和本报告。方法包括人工推导、独立 SymPy 符号核查、复跑原 ODE 检查及原始文献的模型抽查。未使用 Lean 4；这些结果不是 Lean 内核验证过的 PDE 证明。

## 结论与证据对应

| 研究结论 | 评估 | 位置 |
| --- | --- | --- |
| 一维 Euler 不变子类及简单波首次激波时间 | 推导正确；证明的是该不变子类 | [analysis.tex:40](D:/Code/hxh_math/research/analysis.tex:40) |
| 给定二维初值立即产生水平平均声学源 | 正确；该零模并非全空间总动量 | [analysis.tex:130](D:/Code/hxh_math/research/analysis.tex:130) |
| 线性谱及无阻尼声学核 | 正确；水平零模还含剪切分量 | [analysis.tex:103](D:/Code/hxh_math/research/analysis.tex:103) |
| 一阶密度为零和二阶全场公式 | 独立系数核查通过 | [analysis.tex:140](D:/Code/hxh_math/research/analysis.tex:140) |
| 指定初值的二阶残留声学非零 | 恒等式及散射积分核查通过 | [analysis.tex:184](D:/Code/hxh_math/research/analysis.tex:184) |
| 精确平均方程及压力修正 | 结构正确；应力不能任意删除 | [analysis.tex:225](D:/Code/hxh_math/research/analysis.tex:225) |
| 对数时间交接及寿命下界 | 未发现闭合障碍，需补全局部构造和差值引理 | [program.tex:4](D:/Code/hxh_math/research/program.tex:4) |
| 任意固定激波前慢时间上的 Burgers 逼近 | 尚未证明 | [program.tex:153](D:/Code/hxh_math/research/program.tex:153) |
| 真正二维激波形成及寿命渐近 | 尚未证明；有效 Burgers 的爆破不自动传递给原 PDE | [program.tex:186](D:/Code/hxh_math/research/program.tex:186) |

这里“核查通过”仅表示所检查的恒等式成立，不表示整份研究没有其他错误。

## 需要补齐的证明

**第一处是局部适定性和延拓的完整表述，重要程度中等。** [program.tex:68](D:/Code/hxh_math/research/program.tex:68) 用一句话概述正则化、紧性和唯一性；[program.tex:100](D:/Code/hxh_math/research/program.tex:100) 直接给出差值能量估计。作为研究备忘录可以接受，作为自足的定理证明还应给出：解空间及时间连续性、统一正密度区间、正则化的一致估计和取极限、唯一性、延拓准则，以及高阶光滑性传播。未经过同行评审本身不是数学漏洞；这里的问题是必要论证仍被压缩为概述。

**第二处是慢时间主定理尚未闭合，重要程度高。** 寿命下界

\[
T_{\max}(\varepsilon)\ge t_\varepsilon+c_*\varepsilon^{-2}
\]

只给出某个正常数 \(c_*\)，不能替代对每个 \(0<\tau_0<\tau_*\) 的存在性和一致逼近。若沿稿件的平均方程路线证明，至少要处理 [program.tex:176](D:/Code/hxh_math/research/program.tex:176) 所列的应力累计与双声学族交叉项。例如，逐时刻 \(\mathcal S_y=O(\varepsilon^4)\) 在长度 \(O(\varepsilon^{-2})\) 的区间上只能给出 \(O(\varepsilon^2)\) 累计误差，已与主波同阶；需要更小的累计界或有效的抵消论证。跨声学族的零均值观察也需要正常形或平均化误差估计，不能代替它。

**第三处是从激波前逼近到原方程爆破缺少独立论证，重要程度高。** 即使补完上一层定理，对每个固定 \(\tau_0<\tau_*\) 成立的误差界也可能在 \(\tau_0\uparrow\tau_*\) 时失控。必须另行控制压缩梯度或逆特征密度，并处理水平黏性在变形声学坐标中的二阶项。这正是 [program.tex:208](D:/Code/hxh_math/research/program.tex:208) 和 [program.tex:213](D:/Code/hxh_math/research/program.tex:213) 已指出、但尚未解决的问题。

因此，当前可保留“有效 Burgers 的形式激波慢时间为
\(8\sqrt{\mu^2+A\gamma}/(\gamma+1)\)”；不能把它改写为原二维系统已证的首次激波时间。

## 能量论证没有出现所怀疑的导数损失

令 \(X=\|(q,u)\|_{H^m}\)、\(D_m=\|\partial_xu\|_{H^m}\)。在密度固定正区间及小解条件下，扩散系数 \(a(q)=\mu/\rho(q)\) 满足

\[
\|\nabla a\|_\infty+\|a-\mu\|_{H^m}\lesssim X,
\qquad \|\partial_x^2u\|_\infty\lesssim X.
\]

稿件取 \(m\ge6\)，足够支持上述 Sobolev 嵌入。将其扩散交换子与 \(D^\alpha u\) 配对，得到

\[
|\text{交换子能量项}|
\lesssim X^2D_m+X^3
\le \eta D_m^2+C_\eta X^4+CX^3.
\]

主扩散项分部积分产生的 \(a_x\) 项也能这样吸收；对称双曲项由 \(CX^3\) 控制。因此

\[
E_m'+c_\mu D_m^2\lesssim E_m^{3/2}+E_m^2
\]

在这里没有要求不可获得的垂直高阶耗散。不能仅因黏性只有水平方向，就认定 [program.tex:49](D:/Code/hxh_math/research/program.tex:49) 错误。

差值扩散项应明确写成

\[
a(q)\partial_x^2e_u+
\bigl(a(q)-a(q^a)\bigr)\partial_x^2u^a.
\]

第二项的高阶导数由已知近似解承担：
\(\|u^a\|_{H^{m+2}}=O(\varepsilon)\)，故其 \(H^m\) 范数为 \(O(\varepsilon\|e\|_{H^m})\)。这给出补全稿件差值估计的具体办法。更高阶传播可使用系数由低阶 \(X\) 控制的 tame 估计，不必为每个导数阶数重新要求更小的振幅。

在 \(t_\varepsilon=L\log(1/\varepsilon)\)、\(L>1/\mu\) 时，一阶余波除以 \(\varepsilon^2\) 后为 \(\varepsilon^{\mu L-1}\to0\)。二阶水平模式只需趋零，并不需要另加涉及其慢衰减率的 \(L\) 条件。二阶误差除以 \(\varepsilon^2\) 后也趋零，故交接时 \(X=O(\varepsilon^2)\)。从 \(X'\lesssim X^2+X^3\) 重启，确实产生 \(c_*\varepsilon^{-2}\) 的寿命下界。

这些计算支持原证明路线；完整的正则化与延拓引理仍应补写。

## 一个精确抵消例子划定推广边界

“瞬时产生声学源”不保证“留下非零二阶自由声学波”。取 \(\mu=c=1\)，例如令 \(A=1/\gamma\)、\(\gamma>1\)，以及

\[
\begin{aligned}
W_2={}&\sqrt5\cos x\cos y+2\sqrt{17}\cos(2x)\sin y
       +\sqrt{123}\cos(3x)\cos y,\\
W_1={}&\sqrt5\sin x\sin y-\sqrt{17}\sin(2x)\cos y
       +\frac{\sqrt{123}}3\sin(3x)\sin y.
\end{aligned}
\]

它光滑、真正依赖二维坐标、散度为零、水平平均为零。以 \(\rho_0=1\)、\(u_0=\varepsilon W\) 初始化，其二阶平均源为

\[
F(t,y)=\frac12\bigl(5e^{-2t}-68e^{-8t}+123e^{-18t}\bigr)\sin(2y),
\qquad F(0,y)=30\sin(2y)\ne0.
\]

但两支二次散射轮廓都等于零，因为

\[
\frac5{1-i}-\frac{68}{4-i}+\frac{123}{9-i}
=\frac52(1+i)-4(4+i)+\frac32(9+i)=0,
\]

另一族取共轭。新增脚本同时用源项积分和复数乘子核查这一抵消。

这给出了 [analysis.tex:300](D:/Code/hxh_math/research/analysis.tex:300) 所述警告的具体例子。**它不是原稿指定单模初值的反例，也不证明这个新初值的完整 PDE 没有高阶声学波或永不爆破。** 它仅否定把单模非抵消结论无条件推广到一般涡旋的做法。

## 文献适用范围抽查

直接核查原始来源后，稿件将下列结果与本模型区分的做法有依据：

| 来源 | 原文所研究的结构 | 对本项目的意义 |
| --- | --- | --- |
| [Gao 等式 1.1](https://arxiv.org/html/2605.23233v1) | 半空间中包含 \(\Delta_hu+\varepsilon\partial_3^2u+\nabla\operatorname{div}u\) | 令垂直黏性参数为零仍保留散度黏性，不能直接覆盖本模型 |
| [Feng 等式 1.1](https://arxiv.org/html/2607.06264v1) | \(\mu\Delta_hu+(\mu+\zeta)\nabla\operatorname{div}u\)，且 \(3\zeta+2\mu>0\) | 有 \(\mu+\zeta>\mu/3>0\)，其声学耗散机制与本题不同 |
| [Masmoudi 等式 1.1](https://arxiv.org/html/2501.05063v1) | 低 Mach 数、各向异性黏性和 Dirichlet 边界 | 可参考近似解与余项方法，不能直接视为固定声速的长期定理 |
| [Huang 等式 1.1 至 1.2](https://arxiv.org/html/2409.01142v1) | Couette 背景及全拉普拉斯、散度黏性 | 零模含扩散，与本题的无阻尼声学零模不同 |
| [Collot 等式 1.1 和 1.6 至 1.7](https://arxiv.org/html/1803.07826v2) | 标量横向黏性 Burgers；奇对称轴上有闭合方程 | 不能直接替代密度与两个速度分量的耦合估计 |
| [Luk 与 Speck 原始摘要](https://arxiv.org/abs/1610.00737) | 无黏可压缩 Euler 的平面激波扰动，允许小涡量 | 不含本题的方向退化二阶黏性项 |

这次是模型和适用范围抽查，不是对全部参考文献的所有定理条件逐条完成全文复审，也不是穷尽性新颖性检索。未发现这些文献直接推翻当前模型的二阶计算；同样不能由此认证该研究方向具有首创性。

## 可复现的验证结果

[check_algebra.py](D:/Code/hxh_math/research/audit/check_algebra.py) 新增 **58 项**精确符号检查，全部通过；并导入原脚本复跑其 **13 项**符号检查和 **4 组** RK4 检查，全部通过。RK4 最大绝对误差为约 \(7.4083\times10^{-12}\)。

新增覆盖包括一般一二阶保守式与原始式的系数一致性、指定初值的全场系数方程、近似解三四阶残差、三阶平均残差消失、四阶压力修正、水平阻尼分支初值和临界极限、散射符号、黎曼速度系数及上述抵消例子。原脚本只检验水平 \(H\) 的 ODE，没有覆盖其初值及全场残差；本次补上了这些覆盖缺口。

在项目根目录运行：

```powershell
python research/audit/check_algebra.py
```

环境为 Python 3.14.0、SymPy 1.14.0。结果写入 [algebra_results.json](D:/Code/hxh_math/research/audit/algebra_results.json)，其中记录环境版本、每项化简余式、原检查复跑结果以及脚本和研究源文件 SHA256。新增检查失败会返回非零退出码。原 `verification.json` 未改写。

数值部分只积分二阶平均 ODE，不求解原 PDE。符号部分也不验证 Sobolev 不等式、存在性、极限交换或奇性形成。当前环境的 PATH 未找到 `lean`、`lake`、`elan`；本次采用上述方法，未建立 Lean 工程。

## 建议的下一步

优先把局部适定性、差值稳定性和高阶延拓整理为完整引理，使对数交接与寿命下界成为可以独立审查的结论。随后再完成激波前慢时间逼近，最后处理真正二维激波形成。

慢时间部分可考虑在交接时初始化一个精确的一维非线性 Euler 背景，再估计完整二维解与该背景的差值。若背景在固定激波前慢时间内的高阶扰动为 \(O(\varepsilon^2)\)，且差值稳定估计能以该量为增长系数，就可能传播交接误差，而不必先完成所有非零水平模的精细指数衰减。这是待实施的证明路线，不是本次新增的已证定理。
