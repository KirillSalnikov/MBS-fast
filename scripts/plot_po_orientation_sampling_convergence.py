#!/usr/bin/env python3
"""Plot orientation-grid convergence against the Q=2 result."""

from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


MUELLER_NAMES = [f"M{i}{j}" for i in range(1, 5) for j in range(1, 5)]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    return parser.parse_args()


def read_result(run: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    path = run / "result" / "result.dat"
    table = np.loadtxt(path, skiprows=1)
    return table[:, 0], table[:, 1], table[:, 2:18]


def read_run_metadata(run: Path) -> tuple[int, float]:
    stdout = (run / "stdout.log").read_text(errors="replace")
    match = re.search(r"Orientation grid:\s*(\d+)\s*x\s*(\d+)\s*=\s*(\d+)", stdout)
    if not match:
        match = re.search(r"Random grid:\s*(\d+)\s*x\s*(\d+)\s*=\s*(\d+)", stdout)
    orientations = int(match.group(3)) if match else 0
    resources = (run / "resources.txt").read_text(errors="replace")
    elapsed_match = re.search(r"elapsed_seconds=([0-9.eE+-]+)", resources)
    elapsed = float(elapsed_match.group(1)) if elapsed_match else math.nan
    return orientations, elapsed


def weighted_norm(values: np.ndarray, weights: np.ndarray) -> float:
    return float(np.sqrt(np.sum(weights[:, None] * values * values)))


def normalized_mueller(mueller: np.ndarray) -> np.ndarray:
    m11 = mueller[:, :1]
    if np.any(m11 <= 0):
        raise RuntimeError("M11 must be positive for pointwise Mueller normalization")
    return mueller / m11


def main() -> None:
    root = parse_args().root.resolve()
    runs = []
    for run in root.glob("q_*"):
        if not (run / "result" / "result.dat").is_file():
            continue
        q = float(run.name.removeprefix("q_").replace("p", "."))
        runs.append((q, run))
    runs.sort()
    if not runs or runs[-1][0] != 2.0:
        raise SystemExit("A complete Q=2.0 reference result is required")

    theta_ref, weights_ref, mueller_ref = read_result(runs[-1][1])
    weights_ref = weights_ref / np.sum(weights_ref)
    ref_matrix_norm = weighted_norm(mueller_ref, weights_ref)
    ref_m11_norm = weighted_norm(mueller_ref[:, :1], weights_ref)
    normalized_ref = normalized_mueller(mueller_ref)

    rows = []
    element_global_errors = []
    results_by_q = {}
    for q, run in runs:
        theta, weights, mueller = read_result(run)
        if theta.shape != theta_ref.shape or not np.allclose(theta, theta_ref, rtol=0, atol=1e-9):
            raise RuntimeError(f"Scattering-angle grid differs for Q={q}")
        diff = mueller - mueller_ref
        normalized = normalized_mueller(mueller)
        normalized_diff = normalized - normalized_ref
        results_by_q[q] = mueller
        orientations, elapsed = read_run_metadata(run)
        full_error = 100.0 * weighted_norm(diff, weights_ref) / ref_matrix_norm
        m11_error = 100.0 * weighted_norm(diff[:, :1], weights_ref) / ref_m11_norm
        m11_pointwise = 100.0 * float(np.sqrt(np.mean((mueller[:, 0] / mueller_ref[:, 0] - 1.0) ** 2)))
        ratio_rms = 100.0 * float(np.sqrt(np.mean(normalized_diff[:, 1:] ** 2)))
        m12_ratio_rms = 100.0 * float(np.sqrt(np.mean(normalized_diff[:, 1] ** 2)))
        element_errors = 100.0 * np.sqrt(np.mean(normalized_diff * normalized_diff, axis=0))
        element_errors[0] = m11_pointwise
        element_global_errors.append(element_errors)
        row = {
            "Q": q,
            "orientations": orientations,
            "elapsed_seconds": elapsed,
            "M11_L2_percent": m11_error,
            "Mueller_Frobenius_L2_percent": full_error,
            "M11_pointwise_RMS_percent": m11_pointwise,
            "normalized_Mueller_RMS_percentage_points": ratio_rms,
            "M12_over_M11_RMS_percentage_points": m12_ratio_rms,
        }
        row.update({f"{name}_normalized_RMS_percentage_points": value for name, value in zip(MUELLER_NAMES, element_errors)})
        rows.append(row)

    csv_path = root / "orientation_convergence_metrics.csv"
    with csv_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    q_values = np.array([row["Q"] for row in rows])
    orientations = np.array([row["orientations"] for row in rows])
    elapsed = np.array([row["elapsed_seconds"] for row in rows])
    m11_error = np.array([row["M11_L2_percent"] for row in rows])
    full_error = np.array([row["Mueller_Frobenius_L2_percent"] for row in rows])
    m11_pointwise = np.array([row["M11_pointwise_RMS_percent"] for row in rows])
    ratio_rms = np.array([row["normalized_Mueller_RMS_percentage_points"] for row in rows])
    m12_ratio_rms = np.array([row["M12_over_M11_RMS_percentage_points"] for row in rows])
    heatmap = np.asarray(element_global_errors).T

    plt.rcParams.update({"font.size": 10, "axes.grid": True, "grid.alpha": 0.25})
    fig, axes = plt.subplots(1, 3, figsize=(14.4, 4.6), constrained_layout=True)

    ax = axes[0]
    ax.plot(q_values, orientations, "o-", color="#145da0", label="Число ориентаций")
    ax.set_xlabel("Q ориентационной сетки")
    ax.set_ylabel("Число ориентаций")
    ax.set_yscale("log")
    time_ax = ax.twinx()
    time_ax.plot(q_values, elapsed, "s--", color="#d1495b", label="Время")
    time_ax.set_ylabel("Время, с")
    time_ax.set_yscale("log")
    handles = ax.lines + time_ax.lines
    ax.legend(handles, [line.get_label() for line in handles], loc="upper left")

    ax = axes[1]
    positive = q_values < 2.0
    ax.semilogy(q_values[positive], m11_pointwise[positive], "o-", color="#145da0", label="M11, поточечно")
    ax.semilogy(q_values[positive], ratio_rms[positive], "s-", color="#d1495b", label="Все Mij/M11")
    ax.semilogy(q_values[positive], m12_ratio_rms[positive], "^-", color="#2a9d8f", label="M12/M11")
    ax.axhline(1.0, color="#444444", linestyle=":", linewidth=1, label="1 %")
    ax.axhline(5.0, color="#777777", linestyle="--", linewidth=1, label="5 %")
    ax.set_xlabel("Q ориентационной сетки")
    ax.set_ylabel("Взвешенная L2-ошибка, %")
    ax.legend()

    ax = axes[2]
    image = ax.imshow(
        np.maximum(heatmap[:, :-1], 1e-5),
        aspect="auto",
        origin="lower",
        norm="log",
        extent=[q_values[0] - 0.05, q_values[-2] + 0.05, -0.5, 15.5],
        cmap="viridis",
    )
    ax.grid(False)
    ax.set_xlabel("Q ориентационной сетки")
    ax.set_yticks(range(16), MUELLER_NAMES)
    ax.set_ylabel("Элемент матрицы")
    colorbar = fig.colorbar(image, ax=ax, pad=0.02)
    colorbar.set_label("СКО нормированного элемента, процентные пункты")

    fig.suptitle(
        "Сходимость ориентационного усреднения: полый столбик L=20 мкм, μ=0.6622\n"
        "λ=0.532 мкм, n=1.3116, эталон Q=2",
        fontsize=12,
    )
    for suffix in ("png", "pdf"):
        fig.savefig(root / f"orientation_convergence.{suffix}", dpi=220)
    plt.close(fig)

    selected_q = [q for q in (0.1, 0.3, 0.5, 1.0, 1.5, 2.0) if q in results_by_q]
    fig, axes = plt.subplots(4, 4, figsize=(14.4, 11.0), sharex=True, constrained_layout=True)
    colors = plt.cm.viridis(np.linspace(0.05, 0.95, len(selected_q)))
    for index, (ax, name) in enumerate(zip(axes.flat, MUELLER_NAMES)):
        for q, color in zip(selected_q, colors):
            mueller = results_by_q[q]
            values = mueller[:, 0] if index == 0 else mueller[:, index] / mueller[:, 0]
            ax.plot(theta_ref, values, color=color, linewidth=1.1, label=f"Q={q:g}")
        ax.set_title(name)
        if index == 0:
            ax.set_yscale("log")
            ax.set_ylabel("M11")
        elif index % 4 == 0:
            ax.set_ylabel("Mij/M11")
        if index >= 12:
            ax.set_xlabel("Угол рассеяния, град.")
        ax.grid(True, alpha=0.22)
    axes[0, 0].legend(ncol=2, fontsize=8)
    fig.suptitle("Сходимость всех элементов по ориентационной сетке: L=20 мкм, μ=0.6622", fontsize=13)
    for suffix in ("png", "pdf"):
        fig.savefig(root / f"orientation_convergence_all_mueller.{suffix}", dpi=220)
    plt.close(fig)

    print(csv_path)


if __name__ == "__main__":
    main()
