#!/usr/bin/env python3
"""Stable Fieller half-envelope; use with paired independent replica means.

Algebraically identical to the quadratic formula, but expressed in X-rY
coordinates to avoid subtracting nearly equal O(M11**4) quantities.
"""
import numpy as np


def fieller(numerator, denominator, t=2.364624251):
    x_samples = np.asarray(numerator, dtype=np.float64)
    y_samples = np.asarray(denominator, dtype=np.float64)
    if x_samples.shape[0] != y_samples.shape[0] or x_samples.shape[0] < 2:
        raise ValueError('Need paired replica arrays with at least two replicas')
    count = x_samples.shape[0]
    x, y = x_samples.mean(axis=0), y_samples.mean(axis=0)
    vy = y_samples.var(axis=0, ddof=1)/count
    a = y*y-t*t*vy
    with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
        ratio = x/y
        residual = x_samples-ratio*y_samples
        centered = residual-residual.mean(axis=0)
        vr = (centered*centered).sum(axis=0)/(count*(count-1))
        cov = (centered*(y_samples-y)).sum(axis=0)/(count*(count-1))
        width = (np.abs(t*t*cov)+np.sqrt(t**4*cov*cov+a*t*t*vr))/a
    ok = (a > 0) & (y > 0) & np.isfinite(width)
    return np.where(ok, width, np.inf)
