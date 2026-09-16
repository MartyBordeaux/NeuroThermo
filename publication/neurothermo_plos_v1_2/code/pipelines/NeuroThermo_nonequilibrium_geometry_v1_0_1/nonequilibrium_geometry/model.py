from __future__ import annotations

import math
import numpy as np

try:
    from numba import njit
except ImportError:  # Allows validation/tests before the server environment is installed.
    def njit(*args, **kwargs):
        if args and callable(args[0]):
            return args[0]
        return lambda function: function


@njit(cache=True)
def drift_vec(x, y, z, J, b, r, s, kappa, a, c, d, x_R):
    return (
        y - a * x * x * x + b * x * x - z + kappa * J,
        c - d * x * x - y,
        r * (s * (x - x_R) - z),
    )


@njit(cache=True)
def _stationary_core(seed, initial, burn_steps, sample_steps, stride, h, J, theta, constants, noise):
    np.random.seed(seed)
    x, y, z = initial
    b, r, s, kappa = theta
    a, c, d, x_R = constants
    Dx, Dy, Dz = noise
    for _ in range(burn_steps):
        fx, fy, fz = drift_vec(x, y, z, J, b, r, s, kappa, a, c, d, x_R)
        x += fx * h + math.sqrt(2.0 * Dx * h) * np.random.randn()
        y += fy * h + math.sqrt(2.0 * Dy * h) * np.random.randn()
        z += fz * h + math.sqrt(2.0 * Dz * h) * np.random.randn()
        if not (np.isfinite(x) and np.isfinite(y) and np.isfinite(z)) or max(abs(x), abs(y), abs(z)) > 1e6:
            return np.empty((0, 3)), np.array((x, y, z)), False
    saved = sample_steps // stride
    out = np.empty((saved, 3), dtype=np.float64)
    index = 0
    for step in range(sample_steps):
        fx, fy, fz = drift_vec(x, y, z, J, b, r, s, kappa, a, c, d, x_R)
        x += fx * h + math.sqrt(2.0 * Dx * h) * np.random.randn()
        y += fy * h + math.sqrt(2.0 * Dy * h) * np.random.randn()
        z += fz * h + math.sqrt(2.0 * Dz * h) * np.random.randn()
        if not (np.isfinite(x) and np.isfinite(y) and np.isfinite(z)) or max(abs(x), abs(y), abs(z)) > 1e6:
            return out[:index], np.array((x, y, z)), False
        if (step + 1) % stride == 0:
            out[index] = (x, y, z)
            index += 1
    return out, np.array((x, y, z)), True


def stationary_samples(seed: int, theta, J: float, cfg: dict, start=None):
    model, stationary = cfg["model"], cfg["stationary"]
    dt = float(stationary["dt_ms"])
    h = dt / float(model["model_time_scale_ms"])
    initial = np.asarray(start if start is not None else [model["x0"], model["y0"], model["z0"]], float)
    constants = np.asarray([model["a"], model["c"], model["d"], model["x_R"]], float)
    noise = np.asarray(cfg["noise"]["D"], float) * float(cfg["noise"].get("multiplier", 1.0))
    burn_steps = int(round(float(stationary["burn_ms"]) / dt))
    sample_steps = int(round(float(stationary["sample_ms"]) / dt))
    stride = max(1, int(round(float(stationary["sample_stride_ms"]) / dt)))
    return _stationary_core(int(seed), initial, burn_steps, sample_steps, stride, h, float(J), np.asarray(theta, float), constants, noise)


def drift_on_grid(X, Y, Z, J, theta, cfg):
    b, r, s, kappa = map(float, theta)
    model = cfg["model"]
    return (
        Y - float(model["a"]) * X**3 + b * X**2 - Z + kappa * float(J),
        float(model["c"]) - float(model["d"]) * X**2 - Y,
        r * (s * (X - float(model["x_R"])) - Z),
    )


def fixed_points(theta, J, cfg):
    """All real equilibria and Jacobian eigenvalues of the deterministic HR system."""
    b, r, s, kappa = map(float, theta)
    model = cfg["model"]
    a, c, d, x_R = map(float, (model["a"], model["c"], model["d"], model["x_R"]))
    roots = np.roots([-a, b - d, -s, c + s * x_R + kappa * float(J)])
    result = []
    for root in roots:
        if abs(root.imag) > 1e-8:
            continue
        x = float(root.real)
        y, z = c - d * x * x, s * (x - x_R)
        jacobian = np.array([
            [-3.0 * a * x * x + 2.0 * b * x, 1.0, -1.0],
            [-2.0 * d * x, -1.0, 0.0],
            [r * s, 0.0, -r],
        ])
        eig = np.linalg.eigvals(jacobian)
        result.append((x, y, z, eig))
    return result


@njit(cache=True)
def _deterministic_trace(initial, steps, retain, h, J, theta, constants):
    x, y, z = initial
    b, r, s, kappa = theta
    a, c, d, x_R = constants
    output = np.empty((retain, 3))
    start = steps - retain
    out_index = 0
    for index in range(steps):
        fx, fy, fz = drift_vec(x, y, z, J, b, r, s, kappa, a, c, d, x_R)
        x += h * fx
        y += h * fy
        z += h * fz
        if not (np.isfinite(x) and np.isfinite(y) and np.isfinite(z)) or max(abs(x), abs(y), abs(z)) > 1e6:
            return output[:out_index], False
        if index >= start:
            output[out_index] = (x, y, z)
            out_index += 1
    return output, True


def deterministic_trace(theta, J, cfg):
    pcfg, model = cfg["preflight"], cfg["model"]
    dt = float(pcfg["dt_ms"])
    h = dt / float(model["model_time_scale_ms"])
    steps = int(round((float(pcfg["burn_ms"]) + float(pcfg["sample_ms"])) / dt))
    retain = int(round(float(pcfg["sample_ms"]) / dt))
    initial = np.asarray([model["x0"], model["y0"], model["z0"]], float)
    constants = np.asarray([model["a"], model["c"], model["d"], model["x_R"]], float)
    return _deterministic_trace(initial, steps, retain, h, float(J), np.asarray(theta, float), constants)
