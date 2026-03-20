#!/usr/bin/env python3
"""
Generate optimal non-uniform theta grid for MBS-raw scattering computations.

For large particles (x >> 1), the forward scattering peak has angular width
~lambda/D = pi/x radians. A uniform grid with Delta_theta ~ 0.1 deg wastes
most computation on the featureless region theta > 30 deg.

This script creates a three-zone grid:
  1. Forward peak (0 to ~5*lambda/D): fine spacing = lambda/(2*D) = pi/(2*x) rad
  2. Transition (~5*lambda/D to 10 deg): logarithmic spacing
  3. Background (10 to 180 deg): coarse spacing = 2 deg

Usage:
    python3 generate_theta_grid.py <size_parameter_x> [output_file]

    If output_file is not specified, prints to stdout.
    size_parameter_x = pi * D / lambda

Examples:
    python3 generate_theta_grid.py 100 > /tmp/tgrid_x100.txt
    python3 generate_theta_grid.py 500 /tmp/tgrid_x500.txt
"""

import sys
import numpy as np


def generate_theta_grid(x, fine_factor=0.1, n_transition=20, coarse_step=2.0):
    """
    Generate non-uniform theta grid for size parameter x.

    Parameters
    ----------
    x : float
        Size parameter (pi * D / lambda).
    fine_factor : float
        Fine step = fine_factor * (180/x) degrees. Default 0.1 gives ~10 pts per peak.
    n_transition : int
        Number of points in transition zone (log-spaced).
    coarse_step : float
        Step in degrees for the background region.

    Returns
    -------
    theta_deg : numpy array
        Theta values in degrees, sorted, starting from 0.
    """
    # Forward peak width (degrees)
    peak_width_deg = 180.0 / x  # = pi/x in radians, converted to degrees

    # Fine step (degrees)
    fine_step = fine_factor * peak_width_deg
    fine_step = max(fine_step, 0.01)  # minimum 0.01 deg
    fine_step = min(fine_step, 1.0)   # maximum 1 deg (for small x)

    # Zone boundaries
    fine_end = max(5.0 * peak_width_deg, 1.0)  # at least 1 degree
    fine_end = min(fine_end, 10.0)               # at most 10 degrees
    transition_end = 10.0                         # transition zone ends at 10 deg

    if fine_end >= transition_end:
        # No transition zone needed
        transition_end = fine_end

    # Zone 1: Fine grid (0 to fine_end)
    fine_points = np.arange(0, fine_end, fine_step)
    if len(fine_points) == 0 or fine_points[-1] < fine_end:
        fine_points = np.append(fine_points, fine_end)

    # Zone 2: Transition (fine_end to transition_end, log-spaced)
    if transition_end > fine_end + fine_step:
        transition_points = np.geomspace(fine_end, transition_end, n_transition + 1)[1:]
    else:
        transition_points = np.array([])

    # Zone 3: Coarse grid (transition_end to 180)
    coarse_points = np.arange(transition_end + coarse_step, 180.0 + 0.01, coarse_step)
    if len(coarse_points) == 0 or abs(coarse_points[-1] - 180.0) > 0.01:
        coarse_points = np.append(coarse_points, 180.0)

    # Merge all zones
    all_points = np.concatenate([fine_points, transition_points, coarse_points])

    # Remove duplicates and sort
    all_points = np.unique(np.round(all_points, decimals=6))

    # Ensure 0 and 180 are present
    if all_points[0] > 0.001:
        all_points = np.insert(all_points, 0, 0.0)
    if all_points[-1] < 179.99:
        all_points = np.append(all_points, 180.0)

    return all_points


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    x = float(sys.argv[1])

    if x <= 0:
        print("Error: size parameter must be positive", file=sys.stderr)
        sys.exit(1)

    theta = generate_theta_grid(x)

    # Count grid zones for info
    n_fine = np.sum(theta <= 5 * 180.0 / x)
    n_total = len(theta)
    equiv_uniform = int(180.0 / (0.5 * 180.0 / x)) + 1 if x > 10 else n_total

    # Print info to stderr
    print(f"# Size parameter x = {x}", file=sys.stderr)
    print(f"# Forward peak width: {180.0/x:.4f} deg", file=sys.stderr)
    print(f"# Grid points: {n_total} (vs {equiv_uniform} for equivalent uniform)", file=sys.stderr)
    print(f"# Reduction factor: {equiv_uniform / n_total:.1f}x", file=sys.stderr)

    # Output
    out = sys.stdout
    if len(sys.argv) >= 3:
        out = open(sys.argv[2], 'w')

    out.write(f"# Non-uniform theta grid for x={x}\n")
    out.write(f"# {n_total} points\n")
    for t in theta:
        out.write(f"{t:.6f}\n")

    if out is not sys.stdout:
        out.close()


if __name__ == "__main__":
    main()
