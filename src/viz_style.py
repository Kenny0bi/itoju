"""Shared figure style for itoju (IEEE two-column print).

The trait palette was validated with the dataviz validator (light mode, paper white): adjacent and
all-pairs CVD and normal-vision floors pass for the five psychiatric traits.
"""
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
FIG_DIR = ROOT / "figures"
FIG_DIR.mkdir(exist_ok=True)

# IEEE column widths (inches)
SINGLE = 3.5
DOUBLE = 7.16

# Psychiatric traits, one hue each. Worst all-pairs CVD dE 10.4, normal-vision dE 16.3.
TRAIT = {
    "ASD": "#6a48c4",       # violet, the focal trait
    "SCZ": "#3a9ad9",       # sky
    "BIP": "#e2562b",       # palm oil
    "MDD": "#e8a300",       # turmeric: 2.17:1 on white, always direct-label
    "PTSD": "#c2408e",      # hibiscus
    "MDD_EHR": "#e8a300",   # same trait, drawn filled
    "MDD_Clin": "#e8a300",  # same trait, drawn hollow
}
HOLLOW = {"MDD_Clin"}

DIV_NEG = "#2a78d6"
DIV_MID = "#f0efec"
DIV_POS = "#e34948"

SEQ = ["#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7", "#3987e5",
       "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b"]

INK = "#0b0b0b"
INK_2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
SURFACE = "#ffffff"


def apply() -> None:
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans"],
        "font.size": 8,
        "axes.titlesize": 8.5,
        "axes.labelsize": 8,
        "xtick.labelsize": 7.5,
        "ytick.labelsize": 7.5,
        "legend.fontsize": 7.5,
        "axes.edgecolor": MUTED,
        "axes.labelcolor": INK_2,
        "axes.linewidth": 0.6,
        "xtick.color": MUTED,
        "ytick.color": MUTED,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": False,
        "grid.color": GRID,
        "grid.linewidth": 0.5,
        "grid.linestyle": "-",
        "lines.linewidth": 1.4,
        "lines.solid_capstyle": "round",
        "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
        "savefig.dpi": 600,
        "pdf.fonttype": 42,
        "legend.frameon": False,
    })


def cmap_seq():
    return mpl.colors.LinearSegmentedColormap.from_list("itoju_seq", SEQ)


def cmap_div():
    return mpl.colors.LinearSegmentedColormap.from_list("itoju_div", [DIV_NEG, DIV_MID, DIV_POS])


def _is_colorbar(ax) -> bool:
    return hasattr(ax, "_colorbar") or ax.get_label() == "<colorbar>"


def _overlap_is_intentional(owner, other) -> bool:
    """Axes deliberately stacked in the same spot, or an inset placed inside its parent, are not a collision.

    Text-on-text collisions inside an inset are still caught by the text overlap check.
    """
    if owner.get_position().bounds == other.get_position().bounds:
        return True
    return owner in getattr(other, "child_axes", []) or other in getattr(owner, "child_axes", [])


def _rgba_set(artists):
    out = []
    for a in artists:
        cols = []
        for getter in ("get_facecolor", "get_color", "get_edgecolor"):
            if hasattr(a, getter):
                try:
                    c = mpl.colors.to_rgba_array(getattr(a, getter)())
                except (ValueError, TypeError):
                    continue
                cols.extend(c)
        alpha = a.get_alpha()
        for c in cols:
            c = c.copy()
            if alpha is not None:
                c[3] = alpha
            if c[3] > 0:
                out.append(c)
    return out


def check_legend_colors(fig, tol: float = 0.03) -> list[str]:
    """Every legend swatch must match (RGB and alpha) a color actually drawn in its axes."""
    problems = []
    legends = [(ax, ax.get_legend()) for ax in fig.axes if ax.get_legend() is not None]
    legends += [(None, leg) for leg in fig.legends]
    for ax, leg in legends:
        drawn_axes = [ax] if ax is not None else fig.axes
        drawn = []
        for a in drawn_axes:
            drawn += _rgba_set(list(a.patches) + list(a.collections) + list(a.lines))
        drawn = np.array(drawn) if drawn else np.zeros((0, 4))
        for handle, text in zip(leg.legend_handles, leg.get_texts()):
            swatches = [c for c in _rgba_set([handle]) if not np.allclose(c[:3], mpl.colors.to_rgb(SURFACE), atol=tol)]
            if not swatches:
                continue
            if not any(np.any(np.all(np.abs(drawn - s) <= tol, axis=1)) for s in swatches):
                problems.append(f"legend color not in plot: {text.get_text()!r}")
    return problems


def check_tick_labels(fig, rel_tol: float = 0.01) -> list[str]:
    """A numeric tick label must state its tick's value (or its absolute value, for mirrored axes).

    Catches formatters that round 0.25 steps to one decimal, so that -0.25 prints as -0.2.
    Labels with units or suffixes (e.g. '0.5M', '10 kb') are not parsed and are skipped.
    An axis drawn in transformed units declares it with `axis.itoju_value_of = f`, where f maps
    a tick position to the value its label should state (e.g. np.square for sqrt-spaced R^2).
    """
    problems = []
    fig.canvas.draw()
    for ax in fig.axes:
        for axis, lim in ((ax.xaxis, ax.get_xlim()), (ax.yaxis, ax.get_ylim())):
            if axis.get_scale() != "linear":
                continue
            value_of = getattr(axis, "itoju_value_of", None)
            lo, hi = sorted(lim)
            span = abs((value_of(hi) - value_of(lo)) if value_of else (hi - lo)) or 1.0
            for loc, label in zip(axis.get_majorticklocs(), axis.get_majorticklabels()):
                if not (lo <= loc <= hi):
                    continue
                txt = label.get_text().replace("−", "-").strip()
                try:
                    val = float(txt)
                except ValueError:
                    continue
                expected = value_of(loc) if value_of else loc
                if abs(val - expected) > rel_tol * span and abs(val - abs(expected)) > rel_tol * span:
                    problems.append(f"tick label misstates value: {label.get_text()!r} at {loc:g}")
    return problems


MIN_FONT_PT = 7.0  # IEEE: 6 pt is the minimum acceptable, 8 pt preferred; figures are drawn at final print size


def check_font_sizes(fig, min_pt: float = MIN_FONT_PT) -> list[str]:
    """Every visible text element must be at least min_pt at the size the figure is saved."""
    problems = []
    for t in fig.findobj(mpl.text.Text):
        if t.get_visible() and t.get_text().strip() and t.get_alpha() != 0 and t.get_fontsize() < min_pt - 1e-6:
            problems.append(f"text below {min_pt:g} pt ({t.get_fontsize():.1f} pt): {t.get_text()[:40]!r}")
    return problems


def check_layout(fig, min_gap_pt: float = 1.0) -> list[str]:
    """Find text that overlaps other text (within min_gap_pt) or runs outside the figure.

    This catches collisions that are easy to miss by eye at screen size. A figure is only
    accepted when this returns no problems AND the rendered PNG has been inspected.
    """
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    # Tick labels whose tick lies outside the axis view are never drawn; skip them.
    hidden = set()
    for ax in fig.axes:
        for axis, (lo, hi) in ((ax.xaxis, sorted(ax.get_xlim())), (ax.yaxis, sorted(ax.get_ylim()))):
            for tick in axis.get_major_ticks() + axis.get_minor_ticks():
                loc = tick.get_loc()
                if loc is None or not (lo - 1e-9 * max(1, abs(lo)) <= loc <= hi + 1e-9 * max(1, abs(hi))):
                    hidden.update({id(tick.label1), id(tick.label2)})
    texts = [t for t in fig.findobj(mpl.text.Text)
             if id(t) not in hidden and t.get_visible() and t.get_text().strip() and t.get_alpha() != 0]
    boxes = [(t, t.get_window_extent(renderer)) for t in texts]
    # rotated text (radial labels) is tested by its true drawn outline, not its upright bounding box,
    # which for a label at 45 degrees is mostly empty space and would report collisions that are not there
    polys = [_text_polygon(t, renderer) for t in texts]
    rotated = [_is_rotated(t) for t in texts]
    pad = min_gap_pt * fig.dpi / 72
    fb = fig.bbox
    problems = []

    # Text sitting on another panel's plotting area (e.g., a y-label running into the next panel).
    plot_axes = [a for a in fig.axes if a.get_visible()]
    for (t, b), poly, rot in zip(boxes, polys, rotated):
        owner = t.axes
        for a in plot_axes:
            if a is owner or _is_colorbar(a) or (owner is not None and _is_colorbar(owner) and a is getattr(owner, "_colorbar_parent", None)):
                continue
            e = a.get_window_extent(renderer)
            # shrink by 2 px so text merely touching a frame edge is not flagged
            if rot:
                frame = np.array([[e.x0 + 2, e.y0 + 2], [e.x1 - 2, e.y0 + 2], [e.x1 - 2, e.y1 - 2], [e.x0 + 2, e.y1 - 2]])
                hit = _polygons_overlap(poly, frame)
            else:
                hit = b.x0 < e.x1 - 2 and e.x0 + 2 < b.x1 and b.y0 < e.y1 - 2 and e.y0 + 2 < b.y1
            if hit and (owner is None or not _overlap_is_intentional(owner, a)):
                problems.append(f"text on another panel: {t.get_text()!r}")
                break

    problems.extend(check_legend_colors(fig))
    problems.extend(check_tick_labels(fig))
    problems.extend(check_font_sizes(fig))
    for (t, _), poly in zip(boxes, polys):
        if poly[:, 0].min() < fb.x0 - 1 or poly[:, 0].max() > fb.x1 + 1 or poly[:, 1].min() < fb.y0 - 1 or poly[:, 1].max() > fb.y1 + 1:
            problems.append(f"outside figure: {t.get_text()!r}")
    for i in range(len(boxes)):
        a = boxes[i][1]
        for j in range(i + 1, len(boxes)):
            b = boxes[j][1]
            if rotated[i] or rotated[j]:
                hit = _polygons_overlap(polys[i], polys[j], pad)
            else:
                hit = a.x0 - pad < b.x1 and b.x0 - pad < a.x1 and a.y0 - pad < b.y1 and b.y0 - pad < a.y1
            if hit:
                problems.append(f"overlap: {boxes[i][0].get_text()!r} <-> {boxes[j][0].get_text()!r}")
    return problems


def _is_rotated(t) -> bool:
    r = t.get_rotation() % 360
    return min(r, 360 - r) > 1e-6


def _text_polygon(t, renderer) -> np.ndarray:
    """Corners of a text's drawn box in display coordinates, following its rotation."""
    b = t.get_window_extent(renderer)
    if not _is_rotated(t):
        return np.array([[b.x0, b.y0], [b.x1, b.y0], [b.x1, b.y1], [b.x0, b.y1]])
    rot = t.get_rotation()
    t.set_rotation(0)
    b0 = t.get_window_extent(renderer)
    t.set_rotation(rot)
    w, h = b0.width, b0.height
    a = np.deg2rad(rot)
    R = np.array([[np.cos(a), -np.sin(a)], [np.sin(a), np.cos(a)]])
    if t.get_rotation_mode() == "anchor":
        # aligned in the unrotated frame, then turned about the anchor point
        x0 = {"left": 0.0, "center": -w / 2, "right": -w}[t.get_horizontalalignment()]
        y0 = {"bottom": 0.0, "baseline": 0.0, "center": -h / 2, "center_baseline": -h / 2, "top": -h}[t.get_verticalalignment()]
        origin = np.asarray(t.get_transform().transform(t.get_unitless_position()), dtype=float)
    else:
        # rotated first, then its upright bounding box is aligned: the box center is the text center
        x0, y0 = -w / 2, -h / 2
        origin = np.array([(b.x0 + b.x1) / 2, (b.y0 + b.y1) / 2])
    corners = np.array([[x0, y0], [x0 + w, y0], [x0 + w, y0 + h], [x0, y0 + h]])
    return corners @ R.T + origin


def _polygons_overlap(p: np.ndarray, q: np.ndarray, pad: float = 0.0) -> bool:
    """Separating axis test for two convex polygons; closer than pad counts as touching."""
    for poly in (p, q):
        for i in range(len(poly)):
            e = poly[(i + 1) % len(poly)] - poly[i]
            n = np.array([-e[1], e[0]])
            length = np.hypot(*n)
            if length == 0:
                continue
            n = n / length
            pp, qq = p @ n, q @ n
            if pp.max() + pad <= qq.min() or qq.max() + pad <= pp.min():
                return False
    return True


def save(fig, name: str) -> None:
    """Run the checks, then write vector PDF and SVG (paper and poster) and a 600 dpi PNG.

    Vector text and lines scale to any poster size; any rasterized layer (dense hexbins,
    heatmaps) is embedded at 600 dpi, above IEEE's 300 dpi requirement for color images.
    """
    problems = check_layout(fig)
    status = "LAYOUT OK" if not problems else f"LAYOUT PROBLEMS ({len(problems)})"
    print(f"[{name}] {status}")
    for p in problems:
        print(f"   {p}")
    fig.savefig(FIG_DIR / f"{name}.pdf", bbox_inches="tight", dpi=600)
    fig.savefig(FIG_DIR / f"{name}.svg", bbox_inches="tight", dpi=600)
    fig.savefig(FIG_DIR / f"{name}.png", bbox_inches="tight", dpi=600)
    plt.close(fig)
