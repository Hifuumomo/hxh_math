研究主题：仅水平黏性的二维等熵可压缩流体
日期：2026-10-30

阅读 main.pdf。LaTeX 主文件为 main.tex，依次包含 analysis.tex、literature.tex、program.tex、references.tex。
文档明确区分直接恒等式、能量估计证明方案、待证明的慢时间调制和激波猜想。

编译：安装含 ctex 和 Fandol 字体的 TeX Live，在 PowerShell 中执行 ./build.ps1。
公式核对：Python 3 和 SymPy 环境中执行 python verify_expansion.py，结果写入 verification.json。
该脚本只核对展开系数；不模拟原 PDE，不提供完整系统激波形成的数值证据。

本轮主要结果：显式二阶声学选择映射、非抵消、精确未闭合平均方程、
对数时间交接的能量证明方案与 epsilon^{-2} 寿命下界、激波前慢时间主定理目标。
参考文献含原始论文链接、核对的定理编号、算子差异和访问局限。
