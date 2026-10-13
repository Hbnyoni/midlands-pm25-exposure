"""Transparent exposure and transferred-risk calculations, without health records."""
import numpy as np


def who_band(values, cutpoints):
    """Upper-inclusive guideline bands, with -1 reserved for missing values."""
    values = np.asarray(values, dtype=float)
    result = np.searchsorted(np.asarray(cutpoints), values, side='left').astype('int16')
    return np.where(np.isfinite(values), result, -1)


def weighted_mean(values, weights):
    values, weights = np.broadcast_arrays(np.asarray(values, dtype=float), np.asarray(weights, dtype=float))
    valid = np.isfinite(values) & np.isfinite(weights) & (weights > 0)
    denominator = weights[valid].sum()
    return float(np.sum(values[valid] * weights[valid]) / denominator) if denominator > 0 else float('nan')


def transferred_rr(delta_concentration, rr_per_10):
    """Relative contrast for an explicitly specified endpoint and exposure window."""
    if rr_per_10 <= 0:
        raise ValueError('The published risk ratio must be positive.')
    return np.exp(np.log(rr_per_10) * np.asarray(delta_concentration, dtype=float) / 10.0)


def lag01_reduction(daily_concentration, target=15.0):
    """Difference between original and daily-capped lag 0-1 day means.

    Axis zero must be a complete consecutive calendar, with NaN for missing days.
    Missing neighbours stay missing; the first day has no invented preceding value.
    """
    daily = np.asarray(daily_concentration, dtype=float)
    reduction = np.maximum(daily - target, 0)
    delta = np.full(reduction.shape, np.nan)
    delta[1:] = (reduction[1:] + reduction[:-1]) / 2.0
    return delta


def longest_exceedance_run(daily_concentration, target=15.0):
    """Consecutive observed exceedances; missing days break rather than fill runs."""
    daily = np.asarray(daily_concentration, dtype=float)
    run = np.zeros(daily.shape[1:], dtype='int32')
    maximum = run.copy()
    for row in daily:
        run = np.where(np.isfinite(row) & (row > target), run + 1, 0)
        maximum = np.maximum(maximum, run)
    return maximum
