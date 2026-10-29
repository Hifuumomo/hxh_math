"""独立代数审计；仅验证展开恒等式，不证明 PDE 的存在性或激波形成。

运行：python research/audit/check_algebra.py
输出：本目录 algebra_results.json。只有全部检查通过才以退出码 0 结束。
"""

import hashlib
import importlib.util
import json
import platform
import sys
from pathlib import Path

import sympy as sp


def audit():
    """从方程提取系数，并与稿件声称的表达式逐项比较。"""
    t, x, y, eps = sp.symbols("t x y eps", real=True)
    mu, c, delta = sp.symbols("mu c delta", positive=True)
    gamma = sp.symbols("gamma", positive=True)
    p0, p1, p2 = sp.symbols("p0 p1 p2", real=True)
    coords = (x, y)
    results = {}

    def check(name, expr):
        reduced = sp.simplify(sp.trigsimp(sp.expand(expr)))
        results[name] = {"passed": reduced == 0, "residual": str(reduced)}

    def div(vector):
        return sum(sp.diff(vector[i], coords[i]) for i in range(2))

    def adv(vector, scalar):
        return sum(vector[i] * sp.diff(scalar, coords[i]) for i in range(2))

    def coefficient(expr, degree):
        return sp.expand(expr).coeff(eps, degree)

    a = sp.Function("a")(t, x, y)
    b = sp.Function("b")(t, x, y)
    U = tuple(sp.Function(f"U{i}")(t, x, y) for i in range(2))
    Z = tuple(sp.Function(f"Z{i}")(t, x, y) for i in range(2))
    rho = 1 + eps * a + eps**2 * b
    velocity = tuple(eps * U[i] + eps**2 * Z[i] for i in range(2))
    # p0,p1,p2 分别代表 p(1),p'(1),p''(1)，避免先假定特定压力律。
    pressure = p0 + p1 * eps * a + eps**2 * (p1 * b + p2 * a**2 / 2)
    mass = sp.diff(rho, t) + div(tuple(rho * v for v in velocity))
    mass1 = sp.diff(a, t) + div(U)
    mass2 = sp.diff(b, t) + div(Z) + div(tuple(a * v for v in U))
    check("generic_mass_order_1", coefficient(mass, 1) - mass1)
    check("generic_mass_order_2", coefficient(mass, 2) - mass2)

    inverse_rho = 1 - eps * a + eps**2 * (a**2 - b)
    for i, coord in enumerate(coords):
        conservative = (
            sp.diff(rho * velocity[i], t)
            + div(tuple(rho * velocity[i] * v for v in velocity))
            + sp.diff(pressure, coord)
            - mu * sp.diff(velocity[i], x, 2)
        )
        primitive = (
            sp.diff(velocity[i], t)
            + adv(velocity, velocity[i])
            + inverse_rho * sp.diff(pressure, coord)
            - mu * inverse_rho * sp.diff(velocity[i], x, 2)
        )
        primitive1 = sp.diff(U[i], t) + p1 * sp.diff(a, coord) - mu * sp.diff(U[i], x, 2)
        primitive2 = (
            sp.diff(Z[i], t) + adv(U, U[i]) + p1 * sp.diff(b, coord)
            + (p2 - p1) * a * sp.diff(a, coord)
            - mu * sp.diff(Z[i], x, 2) + mu * a * sp.diff(U[i], x, 2)
        )
        check(f"generic_primitive_order_1_{i}", coefficient(primitive, 1) - primitive1)
        check(f"generic_primitive_order_2_{i}", coefficient(primitive, 2) - primitive2)
        check(f"conservative_primitive_identity_order_1_{i}", coefficient(conservative, 1) - primitive1)
        check(f"conservative_primitive_identity_order_2_{i}", coefficient(conservative, 2) - primitive2 - a * primitive1 - U[i] * mass1)

    W = (sp.sin(x) * sp.sin(y), sp.cos(x) * sp.cos(y))
    U = tuple(sp.exp(-mu * t) * v for v in W)
    H, B, P, Q = (sp.Function(name)(t) for name in ("H", "B", "P", "Q"))
    b = H * sp.cos(2 * x) + B * sp.cos(2 * y)
    Z = (P * sp.sin(2 * x), Q * sp.sin(2 * y))
    damping = sp.exp(-2 * mu * t)
    ode = {
        sp.diff(H, t): -2 * P,
        sp.diff(B, t): -2 * Q,
        sp.diff(P, t): -4 * mu * P + 2 * c**2 * H - damping / 2,
        sp.diff(Q, t): 2 * c**2 * B + damping / 2,
    }
    check("full_field_linear_mass", div(U))
    for i in range(2):
        check(f"full_field_linear_momentum_{i}", sp.diff(U[i], t) - mu * sp.diff(U[i], x, 2))
    check("full_field_second_mass", (sp.diff(b, t) + div(Z)).subs(ode))
    for i, coord in enumerate(coords):
        residual2 = sp.diff(Z[i], t) + c**2 * sp.diff(b, coord) - mu * sp.diff(Z[i], x, 2) + adv(U, U[i])
        check(f"full_field_second_momentum_{i}", residual2.subs(ode))
    initial_ut = tuple(-adv(W, v) + mu * sp.diff(v, x, 2) for v in W)
    check("initial_density_second_time_derivative", -div(initial_ut) - sp.cos(2 * x) + sp.cos(2 * y))

    # 直接展开二阶近似解代入原始方程后的残差，核对三、四阶源。
    rho_a = 1 + eps**2 * b
    u_a = tuple(eps * U[i] + eps**2 * Z[i] for i in range(2))
    mass_a = sp.diff(rho_a, t) + div(tuple(rho_a * v for v in u_a))
    cubic_mass = div(tuple(b * v for v in U))
    quartic_mass = div(tuple(b * v for v in Z))
    check("approximate_mass_order_3", coefficient(mass_a, 3) - cubic_mass)
    check("approximate_mass_order_4", coefficient(mass_a, 4) - quartic_mass)
    avg_x = lambda expr: sp.simplify(sp.integrate(sp.expand_trig(expr), (x, 0, 2 * sp.pi)) / (2 * sp.pi))
    check("mean_cubic_mass_vanishes", avg_x(cubic_mass))
    for i, coord in enumerate(coords):
        primitive_a = (
            sp.diff(u_a[i], t) + adv(u_a, u_a[i])
            + c**2 * (1 + (gamma - 2) * eps**2 * b) * eps**2 * sp.diff(b, coord)
            - mu * (1 - eps**2 * b + eps**4 * b**2) * sp.diff(u_a[i], x, 2)
        )
        cubic = adv(U, Z[i]) + adv(Z, U[i]) + mu * b * sp.diff(U[i], x, 2)
        quartic = adv(Z, Z[i]) + c**2 * (gamma - 2) * b * sp.diff(b, coord) + mu * b * sp.diff(Z[i], x, 2)
        check(f"approximate_velocity_order_3_{i}", coefficient(primitive_a, 3) - cubic)
        check(f"approximate_velocity_order_4_{i}", coefficient(primitive_a, 4) - quartic)
        check(f"mean_cubic_velocity_residual_{i}_vanishes", avg_x(cubic))
    check("mean_pressure_order_4", p2 * (avg_x(b**2) - avg_x(b)**2) / 2 - p2 * H**2 / 4)
    check("mean_pressure_order_4_has_no_vertical_force", sp.diff(p2 * (avg_x(b**2) - avg_x(b)**2) / 2, y))
    check("mean_reynolds_order_2", avg_x(U[1]**2) - damping * sp.cos(y)**2 / 2)
    check("mean_horizontal_quadratic_flux", avg_x(U[0] * U[1]))

    h_cases = {
        "underdamped": (damping * (1 - sp.cos(2 * delta * t)) / (4 * delta**2), mu**2 + delta**2),
        "critical": (damping * t**2 / 2, mu**2),
        "overdamped": (damping * (sp.cosh(2 * delta * t) - 1) / (4 * delta**2), mu**2 - delta**2),
    }
    for name, (h, c_squared) in h_cases.items():
        p = -sp.diff(h, t) / 2
        check(f"horizontal_{name}_density_initial", h.subs(t, 0))
        check(f"horizontal_{name}_velocity_initial", p.subs(t, 0))
        check(f"horizontal_{name}_first_order_system", sp.diff(p, t) + 4 * mu * p - 2 * c_squared * h + damping / 2)
    for name in ("underdamped", "overdamped"):
        check(f"horizontal_{name}_critical_limit", sp.limit(h_cases[name][0], delta, 0, dir="+") - h_cases["critical"][0])

    # 用源项的 Laplace 积分独立导出自由声学轨道，并核对散射相位符号。
    s, theta = sp.symbols("s theta", real=True)
    D = mu**2 + c**2
    for sign, label in ((1, "plus"), (-1, "minus")):
        scattering = sp.integrate(sp.exp(-2 * mu * s) * sp.sin(2 * (theta + sign * c * s)) / 2, (s, 0, sp.oo))
        expected = (mu * sp.sin(2 * theta) + sign * c * sp.cos(2 * theta)) / (4 * D)
        check(f"scattering_integral_{label}", scattering - expected)
        # W_2=cos x cos y 在四个 Fourier 模上的系数均为 1/4。
        qhat2 = -sp.I * 2 * (sp.Rational(1, 16) + sp.Rational(1, 16)) / (2 * mu - sign * sp.I * 2 * c)
        expected_hat2 = (sign * c - sp.I * mu) / (8 * D)
        check(f"scattering_fourier_mode_2_{label}", qhat2 - expected_hat2)

    jp, jm = sp.symbols("j_plus j_minus", real=True)
    V = (jp + jm) / 2
    acoustic_c = c + (gamma - 1) * (jp - jm) / 4
    kappa, beta = (gamma + 1) / 4, (3 - gamma) / 4
    check("riemann_characteristic_plus", V + acoustic_c - (c + kappa * jp + beta * jm))
    check("riemann_characteristic_minus", V - acoustic_c - (-c + kappa * jm + beta * jp))
    # sin/cos 导数的振幅平方与形式激波时间；此处不证明原 PDE 爆破。
    derivative_amplitude_squared = (mu / (2 * D))**2 + (c / (2 * D))**2
    check("envelope_derivative_amplitude_squared", derivative_amplitude_squared - 1 / (4 * D))
    shock_time = 1 / (kappa * (1 / (2 * sp.sqrt(D))))
    check("formal_burgers_shock_time", shock_time - 8 * sp.sqrt(D) / (gamma + 1))

    # 负对照：瞬时声源非零，并不推出一般初值的最终二次散射非零。
    # 此数据不是稿件的单个 Taylor--Green 测试模式，不反驳其非抵消结论。
    cancel_w = (
        sp.sqrt(5) * sp.sin(x) * sp.sin(y)
        - sp.sqrt(17) * sp.sin(2 * x) * sp.cos(y)
        + sp.sqrt(123) * sp.sin(3 * x) * sp.sin(y) / 3,
        sp.sqrt(5) * sp.cos(x) * sp.cos(y)
        + 2 * sp.sqrt(17) * sp.cos(2 * x) * sp.sin(y)
        + sp.sqrt(123) * sp.cos(3 * x) * sp.cos(y),
    )
    check("cancellation_control_divergence_free", div(cancel_w))
    for i in range(2):
        check(f"cancellation_control_horizontal_mean_{i}", avg_x(cancel_w[i]))
    initial_force = -sp.diff(avg_x(cancel_w[1] ** 2), y)
    check("cancellation_control_nonzero_initial_force", initial_force - 30 * sp.sin(2 * y))
    heat_velocity = (
        sp.sqrt(5) * sp.exp(-2 * s / 2) * sp.cos(x) * sp.cos(y)
        + 2 * sp.sqrt(17) * sp.exp(-8 * s / 2) * sp.cos(2 * x) * sp.sin(y)
        + sp.sqrt(123) * sp.exp(-18 * s / 2) * sp.cos(3 * x) * sp.cos(y)
    )
    force_amplitude = (5 * sp.exp(-2 * s) - 68 * sp.exp(-8 * s) + 123 * sp.exp(-18 * s)) / 2
    check("cancellation_control_heat_forcing", -sp.diff(avg_x(heat_velocity**2), y) - force_amplitude * sp.sin(2 * y))
    for sign, label in ((1, "plus"), (-1, "minus")):
        multiplier_sum = 5 / (1 - sign * sp.I) - 68 / (4 - sign * sp.I) + 123 / (9 - sign * sp.I)
        check(f"cancellation_control_scattering_multiplier_{label}", multiplier_sum)
        scattering_integral = sp.integrate(force_amplitude * sp.sin(2 * (theta + sign * s)), (s, 0, sp.oo))
        check(f"cancellation_control_scattering_integral_{label}", scattering_integral)
    negative_control = {
        "parameters": {"mu": 1, "c": 1},
        "W1": str(cancel_w[0]),
        "W2": str(cancel_w[1]),
        "F_initial": str(sp.trigsimp(initial_force)),
        "F_at_time_s": str(force_amplitude * sp.sin(2 * y)),
        "Q_plus": "0",
        "Q_minus": "0",
        "interpretation": "验证 analysis.tex 第300行的一般数据相位抵消警告；不反驳稿件指定单模的非抵消结论。",
    }
    return results, negative_control


def main():
    checks, negative_control = audit()
    root = Path(__file__).resolve().parents[1]
    sources = [root / name for name in ("analysis.tex", "program.tex", "verify_expansion.py")]
    # 导入只运行检查函数，不触发原脚本的 __main__，保留原 verification.json。
    original_path = root / "verify_expansion.py"
    spec = importlib.util.spec_from_file_location("original_verify_expansion", original_path)
    original = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(original)
    original_checks = original.symbolic_checks()
    original_rk4 = [original.rk4_check(mu, c) for mu, c in [(1, 1), (0.2, 1), (3, 1), (1, 2)]]
    report = {
        "scope": "SymPy 精确符号恒等式核查；不是 Lean 内核证明、PDE 求解或激波证据。",
        "assumptions": ["mu>0, c>0; 稿件另要求 gamma>1", "过阻尼参数化还需 0<delta<mu", "形式展开系数光滑，水平平均按 2pi 周期归一化"],
        "python_version": sys.version,
        "python_executable": sys.executable,
        "sympy_version": sp.__version__,
        "platform": platform.platform(),
        "audit_script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "source_sha256": {path.name: hashlib.sha256(path.read_bytes()).hexdigest() for path in sources},
        "check_count": len(checks),
        "passed_count": sum(item["passed"] for item in checks.values()),
        "checks": checks,
        "original_script_rerun": {"symbolic_checks": original_checks, "independent_ode_checks": original_rk4},
        "negative_control": negative_control,
        "not_verified": ["局部适定性及延拓", "Sobolev 能量和差值估计", "对数时间误差闭合", "epsilon^-2 时间调制", "原二维 PDE 激波或寿命渐近", "一般无限 Fourier 级数的交换求和/积分及范数连续性"],
    }
    destination = Path(__file__).with_name("algebra_results.json")
    destination.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"check_count": report["check_count"], "passed_count": report["passed_count"], "output": str(destination)}, ensure_ascii=False))
    failed = [name for name, item in checks.items() if not item["passed"]]
    if failed:
        raise AssertionError(failed)


if __name__ == "__main__":
    main()
