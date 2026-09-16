from __future__ import annotations

import numpy as np

from .density import assign_states


def transition_matrix(samples, state_edges, lag_samples, pseudocount):
    states, n_states = assign_states(samples, state_edges)
    counts = np.full((n_states, n_states), float(pseudocount), dtype=float)
    if len(states) > lag_samples:
        np.add.at(counts, (states[:-lag_samples], states[lag_samples:]), 1.0)
    transition = counts / counts.sum(axis=1, keepdims=True)
    occupancy = np.bincount(states, minlength=n_states).astype(float) + float(pseudocount)
    occupancy /= occupancy.sum()
    stationary = occupancy.copy()
    for _ in range(100000):
        updated = stationary @ transition
        if np.max(np.abs(updated - stationary)) < 1e-14:
            stationary = updated
            break
        stationary = updated
    stationary /= stationary.sum()
    return transition, stationary, occupancy, states


def detailed_balance_metrics(transition, pi, cfg, empirical_occupancy=None):
    forward = pi[:, None] * transition
    reverse = forward.T
    flux = forward - reverse
    upper = np.triu_indices_from(flux, k=1)
    denominator = np.sum((forward + reverse)[upper])
    db_violation = float(np.sum(np.abs(flux[upper])) / max(denominator, 1e-300))
    mask = (forward > 0) & (reverse > 0)
    entropy = float(np.sum(forward[mask] * np.log(forward[mask] / reverse[mask])))
    affinities = []
    minimum = float(cfg["markov"]["cycle_min_flux"])
    active = np.where(pi >= float(cfg["markov"]["cycle_min_occupancy"]))[0]
    max_states = int(cfg["markov"]["cycle_max_states"])
    if len(active) > max_states:
        active = active[np.argsort(pi[active])[-max_states:]]
    for ai, i in enumerate(active):
        for aj in range(ai + 1, len(active)):
            j = active[aj]
            for k in active[aj + 1:]:
                terms = [forward[i, j], forward[j, k], forward[k, i], forward[j, i], forward[k, j], forward[i, k]]
                if min(terms) < minimum:
                    continue
                affinities.append(float(np.log((terms[0] * terms[1] * terms[2]) / (terms[3] * terms[4] * terms[5]))))
    absolute = np.abs(affinities)
    return {
        "markov_db_violation": db_violation,
        "markov_entropy_per_lag": entropy,
        "markov_stationarity_residual": float(np.max(np.abs(pi @ transition - pi))),
        "markov_empirical_pi_l1": (
            float(np.sum(np.abs(pi - empirical_occupancy))) if empirical_occupancy is not None else np.nan
        ),
        "cycle_count": int(len(affinities)),
        "cycle_affinity_median_abs": float(np.median(absolute)) if len(absolute) else np.nan,
        "cycle_affinity_max_abs": float(np.max(absolute)) if len(absolute) else np.nan,
    }, affinities


def exact_hatano_sasa(transitions, stationary):
    """Exact <exp(-Y)>=1 check for discrete switches followed by relaxation kernels."""
    weight = stationary[0].copy()
    for index in range(len(transitions) - 1):
        ratio = stationary[index + 1] / stationary[index]
        weight = (weight * ratio) @ transitions[index + 1]
    return float(weight.sum())


def simulate_hatano_sasa(transitions, stationary, trajectories, seed):
    rng = np.random.default_rng(int(seed))
    n_states = len(stationary[0])
    states = rng.choice(n_states, size=int(trajectories), p=stationary[0])
    Y = np.zeros(int(trajectories), dtype=float)
    for index in range(len(transitions) - 1):
        old_pi, new_pi = stationary[index], stationary[index + 1]
        Y += -np.log(new_pi[states]) + np.log(old_pi[states])
        cumulative = np.cumsum(transitions[index + 1], axis=1)
        uniforms = rng.random(len(states))
        states = np.fromiter((np.searchsorted(cumulative[state], value, side="right") for state, value in zip(states, uniforms)), int)
        states = np.minimum(states, n_states - 1)
    exponential = np.exp(-Y)
    ess = float(exponential.sum() ** 2 / np.sum(exponential ** 2))
    return {
        "n_trajectories": int(trajectories),
        "mean_exp_minus_Y": float(exponential.mean()),
        "se_exp_minus_Y": float(exponential.std(ddof=1) / np.sqrt(len(exponential))),
        "median_Y": float(np.median(Y)),
        "ess_fraction": ess / len(exponential),
    }


def simulate_path_probability_ift(transitions, stationary, trajectories, seed):
    """Discrete forward/reversed path-probability ratio; this is not physical Crooks work."""
    rng = np.random.default_rng(int(seed))
    n_states = len(stationary[0])
    states = rng.choice(n_states, size=int(trajectories), p=stationary[0])
    initial = states.copy()
    log_ratio = np.log(stationary[0][states])
    history = [states.copy()]
    for index in range(len(transitions) - 1):
        matrix = transitions[index + 1]
        cumulative = np.cumsum(matrix, axis=1)
        uniforms = rng.random(len(states))
        new_states = np.fromiter((np.searchsorted(cumulative[state], value, side="right") for state, value in zip(states, uniforms)), int)
        new_states = np.minimum(new_states, n_states - 1)
        log_ratio += np.log(matrix[states, new_states])
        states = new_states
        history.append(states.copy())
    log_ratio -= np.log(stationary[-1][states])
    for index in range(len(transitions) - 1):
        reverse_matrix = transitions[len(transitions) - 2 - index]
        right = history[len(transitions) - 1 - index]
        left = history[len(transitions) - 2 - index]
        log_ratio -= np.log(reverse_matrix[right, left])
    exponential = np.exp(-log_ratio)
    return {
        "n_trajectories": int(trajectories),
        "mean_exp_minus_sigma": float(exponential.mean()),
        "se_exp_minus_sigma": float(exponential.std(ddof=1) / np.sqrt(len(exponential))),
        "median_sigma": float(np.median(log_ratio)),
        "initial_state_checksum": int(initial.sum()),
    }
