#!/usr/bin/env python3
"""Generate vector figures used by the MBS-fast manuals."""

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import FancyArrowPatch, Polygon, Circle, Rectangle
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


OUT = Path(__file__).resolve().parent / "figures"
ROOT = Path(__file__).resolve().parents[1]

for font_path in font_manager.findSystemFonts():
    if Path(font_path).name.startswith("cmun"):
        font_manager.fontManager.addfont(font_path)

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["CMU Serif", "Computer Modern Roman", "DejaVu Serif"],
    "mathtext.fontset": "custom",
    "mathtext.rm": "CMU Serif",
    "mathtext.it": "CMU Serif:italic",
    "mathtext.bf": "CMU Serif:bold",
    "pdf.fonttype": 42,
    "axes.unicode_minus": False,
})


def setup(name, w=8.0, h=4.4):
    fig, ax = plt.subplots(figsize=(w, h))
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 5.5)
    return fig, ax, OUT / name


def localized_name(stem, russian):
    suffix = "_ru" if russian else ""
    return f"{stem}{suffix}.pdf"


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


def pipeline(russian=False):
    fig, ax, path = setup(localized_name("manual_pipeline", russian), 9.2, 3.0)
    labels = (
        [
            "грани\nчастицы",
            "трассировка\nлучей",
            "полигоны\nпучков",
            "интеграл\nпо рёбрам",
            "сумма\nДжонса",
            "усреднение\nМюллера",
        ]
        if russian
        else [
            "particle\nfacets",
            "GO ray\ntracing",
            "beam\npolygons",
            "PO edge\nintegral",
            "Jones\nsum",
            "Mueller\naverage",
        ]
    )
    xs = np.linspace(0.4, 8.1, len(labels))
    for x, label in zip(xs, labels):
        box(ax, (x, 2.0), label)
    for x1, x2 in zip(xs[:-1], xs[1:]):
        arrow(ax, (x1 + 1.55, 2.4), (x2, 2.4))
    note = (
        "CPU подготавливает лучи и пучки; CPU или GPU вычисляет дифракцию "
        "на сетке направлений"
        if russian
        else "CPU prepares rays and beams; CPU or GPU evaluates diffraction "
             "on the detector grid"
    )
    ax.text(4.8, 1.15, note, ha="center", fontsize=10)
    save(fig, path)


def aperture_edge(russian=False):
    fig, ax, path = setup(localized_name("manual_aperture_edge", russian), 6.8, 4.6)
    poly = np.array([[2.2, 1.1], [4.8, 1.5], [5.3, 3.4], [3.3, 4.1], [1.5, 2.8]])
    ax.add_patch(Polygon(poly, closed=True, facecolor="#d9ecff", edgecolor="#1f5f99", lw=2))
    ax.text(3.4, 2.6, "апертура A" if russian else "aperture A",
            ha="center", fontsize=11)
    for i, (a, b) in enumerate(zip(poly, np.roll(poly, -1, axis=0))):
        mid = (a + b) / 2
        arrow(ax, tuple(a), tuple(b), color="#1f5f99")
        ax.text(mid[0], mid[1] + 0.2, f"e{i}", fontsize=8, color="#1f5f99")
    arrow(ax, (0.7, 0.7), (1.8, 1.45),
          "пучок" if russian else "beam", "#444444")
    arrow(ax, (5.7, 3.8), (6.5, 4.6),
          "дальняя зона" if russian else "far field", "#b23a48")
    formula = (
        r"$F(q)=\int_A e^{i k q\cdot r}\,dA$ вычисляется как сумма вкладов рёбер"
        if russian
        else r"$F(q)=\int_A e^{i k q\cdot r}\,dA$ is evaluated by summing edge terms"
    )
    ax.text(3.4, 0.45, formula, ha="center", fontsize=11)
    save(fig, path)


def nonconvex_clip(russian=False):
    fig, ax, path = setup(
        localized_name("manual_nonconvex_clipping", russian), 8.0, 4.8)
    subj = np.array([[1.0, 1.0], [5.0, 1.1], [5.1, 4.2], [1.2, 4.0]])
    c1 = np.array([[2.0, 0.8], [3.8, 0.9], [3.9, 4.4], [2.1, 4.3]])
    c2 = np.array([[3.2, 2.2], [5.8, 2.1], [5.6, 3.3], [3.4, 3.5]])
    ax.add_patch(Polygon(subj, closed=True, facecolor="#e6f5dc", edgecolor="#2f7d32", lw=2))
    ax.add_patch(Polygon(c1, closed=True, facecolor="#ffdddd", edgecolor="#b23a48", lw=1.8, alpha=0.8))
    ax.add_patch(Polygon(c2, closed=True, facecolor="#ffdddd", edgecolor="#b23a48", lw=1.8, alpha=0.8))
    ax.text(
        1.2, 4.35,
        "полигон текущего пучка" if russian else "subject beam polygon",
        fontsize=9, color="#2f7d32")
    ax.text(
        3.9, 4.65,
        "проекции блокирующих граней\nна плоскость пучка"
        if russian
        else "blocking facets projected\nonto the subject plane",
        ha="center", fontsize=9, color="#b23a48")
    arrow(ax, (6.6, 1.0), (5.5, 1.8),
          "падающий луч" if russian else "incident ray", "#444444")
    arrow(ax, (6.6, 4.5), (5.7, 3.7),
          "направление проекции" if russian else "clip direction", "#444444")
    result = (
        "видимый полигон = исходный полигон минус объединение проекций"
        if russian
        else "visible polygon = subject minus union(projected blockers)"
    )
    ax.text(3.2, 0.35, result, ha="center", fontsize=11)
    save(fig, path)


def aggregate_clip(russian=False):
    fig, ax, path = setup(
        localized_name("manual_aggregate_clipping", russian), 8.0, 4.4)
    centers = [(2.1, 2.6), (3.3, 2.2), (4.5, 2.9), (5.3, 1.8)]
    for i, c in enumerate(centers):
        hexagon = []
        for k in range(6):
            ang = np.pi / 6 + k * np.pi / 3
            hexagon.append([c[0] + 0.8 * np.cos(ang), c[1] + 0.8 * np.sin(ang)])
        ax.add_patch(Polygon(hexagon, closed=True, facecolor="#e9e2ff",
                             edgecolor="#5b4b9a", lw=1.5, alpha=0.95))
        label = f"часть {i+1}" if russian else f"part {i+1}"
        ax.text(c[0], c[1], label, ha="center", va="center", fontsize=8)
    arrow(ax, (0.6, 3.8), (1.3, 3.4),
          "падающий луч" if russian else "incoming ray", "#333333")
    arrow(ax, (6.8, 0.9), (5.9, 1.3),
          "проверка затенения" if russian else "shadow tests", "#b23a48")
    note = (
        "в агрегате грань может быть закрыта гранями другой части"
        if russian
        else "for aggregates, each facet can be clipped by facets from other parts"
    )
    ax.text(4.0, 0.45, note, ha="center", fontsize=11)
    save(fig, path)


def avx512_pack(russian=False):
    fig, ax, path = setup(
        localized_name("manual_avx512_packing", russian), 9.0, 4.8)
    for lane in range(8):
        y = 3.9 - lane * 0.38
        ax.add_patch(Rectangle((0.8, y), 1.0, 0.28, facecolor="#eef4ff", edgecolor="#3b5f9b"))
        label = f"угол {lane}" if russian else f"theta {lane}"
        ax.text(1.3, y + 0.14, label, ha="center", va="center", fontsize=7)
    box(ax, (2.4, 2.4), "загрузка 8\nфаз" if russian else "load 8\nphases", w=1.1)
    box(ax, (4.0, 2.4), "AVX-512\nsin/cos" if russian else "AVX-512\nsincos", w=1.25)
    box(ax, (5.8, 2.4), "вклады\nрёбер" if russian else "edge\nterms", w=1.1)
    box(ax, (7.3, 2.4), "обновление\nДжонса" if russian else "Jones\nupdates", w=1.15)
    arrow(ax, (1.85, 3.0), (2.4, 2.8))
    arrow(ax, (3.5, 2.8), (4.0, 2.8))
    arrow(ax, (5.25, 2.8), (5.8, 2.8))
    arrow(ax, (6.9, 2.8), (7.3, 2.8))
    first_note = (
        "один регистр __m512d содержит 8 чисел двойной точности; AVX2 — 4"
        if russian
        else "one __m512d register holds 8 double lanes; AVX2 fallback uses 4 lanes"
    )
    second_note = (
        "фазы и скалярные произведения группируются по направлениям или размерам"
        if russian
        else "the code batches trigonometric phases and dot products across detector directions or sizes"
    )
    ax.text(4.6, 1.45, first_note, ha="center", fontsize=11)
    ax.text(4.6, 0.95, second_note, ha="center", fontsize=10)
    save(fig, path)


def gpu_atomics(russian=False):
    fig, ax, path = setup(
        localized_name("manual_gpu_atomics", russian), 8.8, 4.6)
    for i in range(5):
        Circle((1.0, 0.9 + i * 0.65), 0.18, facecolor="#d9ecff", edgecolor="#1f5f99").set_clip_on(False)
        ax.add_patch(Circle((1.0, 0.9 + i * 0.65), 0.18, facecolor="#d9ecff", edgecolor="#1f5f99"))
        arrow(ax, (1.2, 0.9 + i * 0.65), (4.0, 2.2), color="#1f5f99", rad=(i - 2) * 0.08)
    ax.add_patch(Rectangle((4.0, 1.55), 1.5, 1.25, facecolor="#fff1cc", edgecolor="#c28a00", lw=2))
    ax.text(4.75, 2.18, "ячейка\nДжонса" if russian else "Jones\ncell",
            ha="center", va="center", fontsize=10)
    box(ax, (6.2, 1.75),
        "матрица\nМюллера" if russian else "Mueller\nconversion",
        w=1.35)
    arrow(ax, (5.5, 2.18), (6.2, 2.18))
    upper = (
        "атомарный путь: потоки пучков складываются в общей ячейке Джонса"
        if russian
        else "atomic path: many beam threads add into one coherent Jones cell"
    )
    lower = (
        "совмещённый путь: один поток накапливает локально и записывает один раз"
        if russian
        else "no-atomic/fused path: one owner thread accumulates locally, then writes once"
    )
    ax.text(2.5, 3.95, upper, ha="center", fontsize=10)
    ax.text(5.0, 0.75, lower, ha="center", fontsize=10)
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


def particle_gallery(russian=False):
    specs = (
        [
            ("hexagonal_column", "1  Столбик", "L=2, D=1"),
            ("bullet", "2  Пуля", "L=2, D=1; вершина автоматически"),
            ("bullet_rosette", "3  Розетка пуль", "L=2, D=1"),
            ("droxtal", "4  Дрокстал", "масштаб 1"),
            ("concave_hexagonal", "10  Полый столбик", "L=2, D=1; впадина 30°"),
            ("two_column_aggregate", "12  Агрегат столбиков", "L=2, D=1; 2 части"),
            ("fixed_aggregate", "999  Заданный агрегат", "масштаб 1"),
        ]
        if russian
        else [
            ("hexagonal_column", "1  Hexagonal column", "L=2, D=1"),
            ("bullet", "2  Bullet", "L=2, D=1; cap=auto"),
            ("bullet_rosette", "3  Bullet rosette", "L=2, D=1; cap=auto"),
            ("droxtal", "4  Droxtal", "scale=1"),
            ("concave_hexagonal", "10  Concave column", "L=2, D=1; cavity=30 deg"),
            ("two_column_aggregate", "12  Two-column aggregate", "L=2, D=1; parts=2"),
            ("fixed_aggregate", "999  Fixed aggregate", "scale=1"),
        ]
    )
    colors = ["#6baed6", "#74c476", "#fd8d3c", "#9e9ac8",
              "#e6550d", "#31a354", "#756bb1"]
    fig = plt.figure(figsize=(10.8, 8.0))
    for index, ((stem, title, params), color) in enumerate(zip(specs, colors), 1):
        particle = read_particle(ROOT / "examples" / "particles" / f"{stem}.particle")
        ax = fig.add_subplot(2, 4, index, projection="3d")
        draw_particle(ax, particle, color)
        kind = (
            ("невыпуклая" if particle["concave"] else "выпуклая")
            if russian
            else ("nonconvex" if particle["concave"] else "convex")
        )
        if particle["aggregate"]:
            kind += ", агрегат" if russian else ", aggregate"
        ax.set_title(f"{title}\n{params}\n{kind}",
                     fontsize=14 if russian else 9, pad=0)

    note = fig.add_subplot(2, 4, 8)
    note.axis("off")
    note_text = (
        "Построено по файлам программы,\n"
        "записанным параметром --save-geometry.\n\n"
        "D — диаметр шестиугольника\nмежду вершинами\n"
        "L — длина призмы или ветви\n"
        "масштаб — равномерное изменение размеров"
        if russian
        else "Rendered from the native files\nwritten by --save-geometry.\n\n"
             "D: vertex-to-vertex hexagon diameter\n"
             "L: prism or branch length\n"
             "scale: uniform geometry scale"
    )
    note.text(
        0.03, 0.72, note_text, va="top",
        fontsize=14 if russian else 10, linespacing=1.3)
    fig.subplots_adjust(left=0.02, right=0.98, bottom=0.03, top=0.96,
                        wspace=0.02, hspace=0.12)
    fig.savefig(
        OUT / localized_name("manual_particle_gallery", russian),
        bbox_inches="tight")
    plt.close(fig)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    pipeline(russian=True)
    aperture_edge(russian=True)
    nonconvex_clip(russian=True)
    aggregate_clip(russian=True)
    avx512_pack(russian=True)
    gpu_atomics(russian=True)
    particle_gallery(russian=False)
    particle_gallery(russian=True)


if __name__ == "__main__":
    main()
