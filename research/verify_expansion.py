"""核对研究备忘录的展开系数；这些检查不模拟原 PDE，也不证明激波。"""

import json
import math
from pathlib import Path

import sympy as sp


def symbolic_checks():
    """通过直接微分验证强迫项、ODE、初值与声学非抵消。"""
    t, x, y = sp.symbols("t x y", real=True)
    mu, c = sp.symbols("mu c", positive=True)
    damping = sp.exp(-2 * mu * t)
    den = 4 * (mu**2 + c**2)
    b = (sp.cos(2 * c * t) - mu / c * sp.sin(2 * c * t) - damping) / den
    q = (mu * sp.cos(2 * c * t) + c * sp.sin(2 * c * t) - mu * damping) / den
    w = (sp.sin(x) * sp.sin(y), sp.cos(x) * sp.cos(y))
    force = tuple(-w[0] * sp.diff(v, x) - w[1] * sp.diff(v, y) for v in w)
    free_b, free_q = b + damping / den, q + mu * damping / den
    identities = {
        "initial_divergence": sp.diff(w[0], x) + sp.diff(w[1], y),
        "horizontal_force": force[0] + sp.sin(2 * x) / 2,
        "vertical_force": force[1] - sp.sin(2 * y) / 2,
        "mass_ode": sp.diff(b, t) + 2 * q,
        "momentum_ode": sp.diff(q, t) - 2 * c**2 * b - damping / 2,
        "initial_density": b.subs(t, 0),
        "initial_velocity": q.subs(t, 0),
        "noncancellation": c**2 * free_b**2 + free_q**2 - 1 / (16 * (mu**2 + c**2)),
        "right_wave": free_q * sp.sin(2 * y) + c * free_b * sp.cos(2 * y)
        - (mu * sp.sin(2 * (y - c * t)) + c * sp.cos(2 * (y - c * t))) / den,
        "left_wave": free_q * sp.sin(2 * y) - c * free_b * sp.cos(2 * y)
        - (mu * sp.sin(2 * (y + c * t)) - c * sp.cos(2 * (y + c * t))) / den,
    }
    delta = sp.symbols("delta", positive=True)
    h_cases = {
        "underdamped": (sp.exp(-2 * mu * t) * (1 - sp.cos(2 * delta * t)) / (4 * delta**2), mu**2 + delta**2),
        "critical": (sp.exp(-2 * mu * t) * t**2 / 2, mu**2),
        "overdamped": (sp.exp(-2 * mu * t) * (sp.cosh(2 * delta * t) - 1) / (4 * delta**2), mu**2 - delta**2),
    }
    for name, (h, c_squared) in h_cases.items():
        identities[name] = sp.diff(h, t, 2) + 4 * mu * sp.diff(h, t) + 4 * c_squared * h - damping
    checked = {name: sp.simplify(sp.trigsimp(expr)) == 0 for name, expr in identities.items()}
    assert all(checked.values()), checked
    return checked


def rk4_check(mu, c, dt=0.001, end=8.0):
    """独立积分二阶平均 ODE，和闭式解逐步比较。"""
    def rhs(time, state):
        b, q = state
        return (-2 * q, 2 * c * c * b + 0.5 * math.exp(-2 * mu * time))

    def shifted(state, slope, step):
        return tuple(a + step * b for a, b in zip(state, slope))

    state, worst = (0.0, 0.0), 0.0
    for index in range(round(end / dt)):
        time = index * dt
        k1 = rhs(time, state)
        k2 = rhs(time + dt / 2, shifted(state, k1, dt / 2))
        k3 = rhs(time + dt / 2, shifted(state, k2, dt / 2))
        k4 = rhs(time + dt, shifted(state, k3, dt))
        state = tuple(a + dt * (b + 2 * d + 2 * e + f) / 6 for a, b, d, e, f in zip(state, k1, k2, k3, k4))
        now, den = (index + 1) * dt, 4 * (mu * mu + c * c)
        exact = (
            (math.cos(2 * c * now) - mu / c * math.sin(2 * c * now) - math.exp(-2 * mu * now)) / den,
            (mu * math.cos(2 * c * now) + c * math.sin(2 * c * now) - mu * math.exp(-2 * mu * now)) / den,
        )
        worst = max(worst, *(abs(a - b) for a, b in zip(state, exact)))
    assert worst < 1e-8, (mu, c, worst)
    return {"mu": mu, "c": c, "step": dt, "end_time": end, "max_absolute_error": worst}


if __name__ == "__main__":
    result = {
        "scope": "二阶系数代数及ODE核对，不是原PDE的激波数值证据",
        "symbolic_checks": symbolic_checks(),
        "independent_ode_checks": [rk4_check(mu, c) for mu, c in [(1, 1), (0.2, 1), (3, 1), (1, 2)]],
    }
    destination = Path(__file__).with_name("verification.json")
    destination.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False, indent=2))
