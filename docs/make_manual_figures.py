#!/usr/bin/env python3
"""Generate vector figures used by the MBS-fast manuals."""

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Polygon, Circle, Rectangle
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


OUT = Path(__file__).resolve().parent / "figures"
ROOT = Path(__file__).resolve().parents[1]


def setup(name, w=8.0, h=4.4):
    fig, ax = plt.subplots(figsize=(w, h))
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 5.5)
    return fig, ax, OUT / name


def arrow(ax, a, b, text=None, color="#333333", rad=0.0):
    p = FancyArrowPatch(a, b, arrowstyle="-|>", mutation_scale=14,
                        lw=1.8, color=color,
                        connectionstyle=f"arc3,rad={rad}")
    ax.add_patch(p)
    if text:
        x = (a[0] + b[0]) * 0.5
        y = (a[1] + b[1]) * 0.5
        ax.text(x, y + 0.18, text, ha="center", va="bottom",
                fontsize=9, color=color)


def box(ax, xy, text, fc="#eef4ff", ec="#3b5f9b", w=1.55, h=0.8):
    r = Rectangle(xy, w, h, facecolor=fc, edgecolor=ec, lw=1.5)
    ax.add_patch(r)
    ax.text(xy[0] + w / 2, xy[1] + h / 2, text, ha="center",
            va="center", fontsize=9)
    return r


def save(fig, path):
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def pipeline():
    fig, ax, path = setup("manual_pipeline.pdf", 9.2, 3.0)
    labels = [
        "particle\nfacets",
        "GO ray\ntracing",
        "beam\npolygons",
        "PO edge\nintegral",
        "Jones\nsum",
        "Mueller\naverage",
    ]
    xs = np.linspace(0.4, 8.1, len(labels))
    for x, label in zip(xs, labels):
        box(ax, (x, 2.0), label)
    for x1, x2 in zip(xs[:-1], xs[1:]):
        arrow(ax, (x1 + 1.55, 2.4), (x2, 2.4))
    ax.text(4.8, 1.15, "CPU prepares rays and beams; CPU or GPU evaluates diffraction on the detector grid",
            ha="center", fontsize=10)
    save(fig, path)


def aperture_edge():
    fig, ax, path = setup("manual_aperture_edge.pdf", 6.8, 4.6)
    poly = np.array([[2.2, 1.1], [4.8, 1.5], [5.3, 3.4], [3.3, 4.1], [1.5, 2.8]])
    ax.add_patch(Polygon(poly, closed=True, facecolor="#d9ecff", edgecolor="#1f5f99", lw=2))
    ax.text(3.4, 2.6, "aperture A", ha="center", fontsize=11)
    for i, (a, b) in enumerate(zip(poly, np.roll(poly, -1, axis=0))):
        mid = (a + b) / 2
        arrow(ax, tuple(a), tuple(b), color="#1f5f99")
        ax.text(mid[0], mid[1] + 0.2, f"e{i}", fontsize=8, color="#1f5f99")
    arrow(ax, (0.7, 0.7), (1.8, 1.45), "beam", "#444444")
    arrow(ax, (5.7, 3.8), (6.5, 4.6), "far field", "#b23a48")
    ax.text(3.4, 0.45, r"$F(q)=\int_A e^{i k q\cdot r}\,dA$ is evaluated by summing edge terms",
            ha="center", fontsize=11)
    save(fig, path)


def nonconvex_clip():
    fig, ax, path = setup("manual_nonconvex_clipping.pdf", 8.0, 4.8)
    subj = np.array([[1.0, 1.0], [5.0, 1.1], [5.1, 4.2], [1.2, 4.0]])
    c1 = np.array([[2.0, 0.8], [3.8, 0.9], [3.9, 4.4], [2.1, 4.3]])
    c2 = np.array([[3.2, 2.2], [5.8, 2.1], [5.6, 3.3], [3.4, 3.5]])
    ax.add_patch(Polygon(subj, closed=True, facecolor="#e6f5dc", edgecolor="#2f7d32", lw=2))
    ax.add_patch(Polygon(c1, closed=True, facecolor="#ffdddd", edgecolor="#b23a48", lw=1.8, alpha=0.8))
    ax.add_patch(Polygon(c2, closed=True, facecolor="#ffdddd", edgecolor="#b23a48", lw=1.8, alpha=0.8))
    ax.text(1.2, 4.35, "subject beam polygon", fontsize=9, color="#2f7d32")
    ax.text(3.9, 4.65, "blocking facets projected\nonto the subject plane", ha="center", fontsize=9, color="#b23a48")
    arrow(ax, (6.6, 1.0), (5.5, 1.8), "incident ray", "#444444")
    arrow(ax, (6.6, 4.5), (5.7, 3.7), "clip direction", "#444444")
    ax.text(3.2, 0.35, "visible polygon = subject minus union(projected blockers)",
            ha="center", fontsize=11)
    save(fig, path)


def aggregate_clip():
    fig, ax, path = setup("manual_aggregate_clipping.pdf", 8.0, 4.4)
    centers = [(2.1, 2.6), (3.3, 2.2), (4.5, 2.9), (5.3, 1.8)]
    for i, c in enumerate(centers):
        hexagon = []
        for k in range(6):
            ang = np.pi / 6 + k * np.pi / 3
            hexagon.append([c[0] + 0.8 * np.cos(ang), c[1] + 0.8 * np.sin(ang)])
        ax.add_patch(Polygon(hexagon, closed=True, facecolor="#e9e2ff",
                             edgecolor="#5b4b9a", lw=1.5, alpha=0.95))
        ax.text(c[0], c[1], f"part {i+1}", ha="center", va="center", fontsize=8)
    arrow(ax, (0.6, 3.8), (1.3, 3.4), "incoming ray", "#333333")
    arrow(ax, (6.8, 0.9), (5.9, 1.3), "shadow tests", "#b23a48")
    ax.text(4.0, 0.45, "for aggregates, each facet can be clipped by facets from other parts",
            ha="center", fontsize=11)
    save(fig, path)


def avx512_pack():
    fig, ax, path = setup("manual_avx512_packing.pdf", 9.0, 4.8)
    for lane in range(8):
        y = 3.9 - lane * 0.38
        ax.add_patch(Rectangle((0.8, y), 1.0, 0.28, facecolor="#eef4ff", edgecolor="#3b5f9b"))
        ax.text(1.3, y + 0.14, f"theta {lane}", ha="center", va="center", fontsize=7)
    box(ax, (2.4, 2.4), "load 8\nphases", w=1.1)
    box(ax, (4.0, 2.4), "AVX-512\nsincos", w=1.25)
    box(ax, (5.8, 2.4), "edge\nterms", w=1.1)
    box(ax, (7.3, 2.4), "Jones\nupdates", w=1.15)
    arrow(ax, (1.85, 3.0), (2.4, 2.8))
    arrow(ax, (3.5, 2.8), (4.0, 2.8))
    arrow(ax, (5.25, 2.8), (5.8, 2.8))
    arrow(ax, (6.9, 2.8), (7.3, 2.8))
    ax.text(4.6, 1.45, "one __m512d register holds 8 double lanes; AVX2 fallback uses 4 lanes",
            ha="center", fontsize=11)
    ax.text(4.6, 0.95, "the code batches trigonometric phases and dot products across detector directions or sizes",
            ha="center", fontsize=10)
    save(fig, path)


def gpu_atomics():
    fig, ax, path = setup("manual_gpu_atomics.pdf", 8.8, 4.6)
    for i in range(5):
        Circle((1.0, 0.9 + i * 0.65), 0.18, facecolor="#d9ecff", edgecolor="#1f5f99").set_clip_on(False)
        ax.add_patch(Circle((1.0, 0.9 + i * 0.65), 0.18, facecolor="#d9ecff", edgecolor="#1f5f99"))
        arrow(ax, (1.2, 0.9 + i * 0.65), (4.0, 2.2), color="#1f5f99", rad=(i - 2) * 0.08)
    ax.add_patch(Rectangle((4.0, 1.55), 1.5, 1.25, facecolor="#fff1cc", edgecolor="#c28a00", lw=2))
    ax.text(4.75, 2.18, "Jones\ncell", ha="center", va="center", fontsize=10)
    box(ax, (6.2, 1.75), "Mueller\nconversion", w=1.35)
    arrow(ax, (5.5, 2.18), (6.2, 2.18))
    ax.text(2.5, 3.95, "atomic path: many beam threads add into one coherent Jones cell",
            ha="center", fontsize=10)
    ax.text(5.0, 0.75, "no-atomic/fused path: one owner thread accumulates locally, then writes once",
            ha="center", fontsize=10)
    save(fig, path)


def read_particle(path):
    """Read the native particle format without inventing geometry for figures."""
    lines = path.read_text(encoding="ascii").splitlines()
    records = []
    cursor = 0
    while len(records) < 3 and cursor < len(lines):
        text = lines[cursor].split("#", 1)[0].strip()
        cursor += 1
        if text:
            records.append(text)
    if len(records) != 3:
        raise ValueError(f"{path}: missing native particle headers")

    facets = []
    current = []
    for line in lines[cursor:]:
        text = line.split("#", 1)[0].strip()
        if not text:
            if current:
                facets.append(np.asarray(current, dtype=float))
                current = []
            continue
        values = [float(value) for value in text.split()]
        if len(values) != 3:
            raise ValueError(f"{path}: vertex must contain x y z")
        current.append(values)
    if current:
        facets.append(np.asarray(current, dtype=float))
    if not facets:
        raise ValueError(f"{path}: no facets")
    return {
        "concave": bool(int(records[0].split()[0])),
        "aggregate": bool(int(records[1].split()[0])),
        "symmetry": tuple(float(value) for value in records[2].split()),
        "facets": facets,
    }


def draw_particle(ax, particle, color):
    points = np.vstack(particle["facets"])
    center = 0.5 * (points.min(axis=0) + points.max(axis=0))
    span = max(np.ptp(points, axis=0).max(), 1e-12)
    facets = [(facet - center) / span for facet in particle["facets"]]
    collection = Poly3DCollection(
        facets, facecolor=color, edgecolor="#20252b", linewidth=0.45,
        alpha=0.78)
    ax.add_collection3d(collection)
    ax.set_xlim(-0.58, 0.58)
    ax.set_ylim(-0.58, 0.58)
    ax.set_zlim(-0.58, 0.58)
    ax.set_box_aspect((1, 1, 1))
    ax.view_init(elev=22, azim=-52)
    ax.set_axis_off()


def particle_gallery():
    specs = [
        ("hexagonal_column", "1  Hexagonal column", "L=2, D=1"),
        ("bullet", "2  Bullet", "L=2, D=1; cap=auto"),
        ("bullet_rosette", "3  Bullet rosette", "L=2, D=1; cap=auto"),
        ("droxtal", "4  Droxtal", "scale=1"),
        ("concave_hexagonal", "10  Concave column", "L=2, D=1; cavity=30 deg"),
        ("two_column_aggregate", "12  Two-column aggregate", "L=2, D=1; parts=2"),
        ("fixed_aggregate", "999  Fixed aggregate", "scale=1"),
    ]
    colors = ["#6baed6", "#74c476", "#fd8d3c", "#9e9ac8",
              "#e6550d", "#31a354", "#756bb1"]
    fig = plt.figure(figsize=(10.8, 8.0))
    for index, ((stem, title, params), color) in enumerate(zip(specs, colors), 1):
        particle = read_particle(ROOT / "examples" / "particles" / f"{stem}.particle")
        ax = fig.add_subplot(2, 4, index, projection="3d")
        draw_particle(ax, particle, color)
        kind = "nonconvex" if particle["concave"] else "convex"
        if particle["aggregate"]:
            kind += ", aggregate"
        ax.set_title(f"{title}\n{params}\n{kind}", fontsize=9, pad=0)

        # Also provide an inspectable raster for each native example.
        single = plt.figure(figsize=(4.0, 4.0))
        single_ax = single.add_subplot(111, projection="3d")
        draw_particle(single_ax, particle, color)
        single_ax.set_title(f"{title}\n{params}", fontsize=11)
        single.savefig(OUT / f"particle_{stem}.png", dpi=220,
                       bbox_inches="tight", facecolor="white")
        plt.close(single)

    note = fig.add_subplot(2, 4, 8)
    note.axis("off")
    note.text(
        0.03, 0.72,
        "Rendered from the native files\nwritten by --save-geometry.\n\n"
        "D: vertex-to-vertex hexagon diameter\n"
        "L: prism or branch length\n"
        "scale: uniform geometry scale",
        va="top", fontsize=10, linespacing=1.45)
    fig.subplots_adjust(left=0.02, right=0.98, bottom=0.03, top=0.96,
                        wspace=0.02, hspace=0.12)
    fig.savefig(OUT / "manual_particle_gallery.pdf", bbox_inches="tight")
    plt.close(fig)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    pipeline()
    aperture_edge()
    nonconvex_clip()
    aggregate_clip()
    avx512_pack()
    gpu_atomics()
    particle_gallery()


if __name__ == "__main__":
    main()
