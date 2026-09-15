"""itoju's figure system: hand-built SVG in print units, with checks that read the rendered PDF.

Every figure is drawn in points at its final printed size (a single IEEE column is 252 pt wide, the
full page 516 pt), so a 7 pt label in the SVG is a 7 pt label in the paper, and the same file is
rasterised at 600 dpi for the README.

The look follows the house style of the other projects (a tracked-caps kicker, a serif headline that
states the finding, monospace labels, a reading panel, a key, a source line) on a ground of its own:
a cool clinical linen, since the subject is the health record. Colour is kept for the psychiatric
traits, validated all-pairs for colour-vision deficiency on this ground (dataviz validator: worst CVD
dE 10.4, normal-vision dE 16.3; sky and turmeric sit below 3:1 contrast, so they always carry a key
entry or a direct label). Gray-green is reserved for noise and reference, never for a result.

Checks, run on the PDF that cairosvg renders (so they see what a reader sees):
  1. no text overlaps other text
  2. no text runs off the canvas
  3. no text below 7 pt
  4. no text sits on a data mark, unless that label was placed on the mark on purpose
  5. every figure declares a key, or states in code why it needs none
"""
import html
import math
from pathlib import Path

import cairosvg
import pymupdf

ROOT = Path(__file__).resolve().parents[1]
FIG_DIR = ROOT / "figures"
FIG_DIR.mkdir(exist_ok=True)

SINGLE = 252.0
DOUBLE = 516.0

GROUND = "#EEF3F1"
PANEL = "#E2EAE7"
PANEL_2 = "#D5DFDB"
GRID = "#CAD6D1"
RULE = "#A9B8B2"
INK = "#14232B"
INK_2 = "#3A4A51"
DIM = "#5E6E72"
NOISE = "#CDD8D4"        # the sleeve: plus or minus 1.96 SE of a difference
NOISE_2 = "#DCE5E1"      # the detectable shift at 80% power, one step lighter

TRAIT = {"ASD": "#6a48c4", "SCZ": "#3a9ad9", "BIP": "#e2562b", "MDD": "#e8a300", "PTSD": "#c2408e",
         "MDD_EHR": "#e8a300", "MDD_Clin": "#e8a300"}
TRAIT_NAME = {"ASD": "autism", "SCZ": "schizophrenia", "BIP": "bipolar", "MDD": "depression", "PTSD": "PTSD",
              "MDD_EHR": "depression, EHR", "MDD_Clin": "depression, clinical"}
FIVE = ["ASD", "SCZ", "BIP", "MDD", "PTSD"]
FIRST = "#c9761a"        # first definition of a pair (validated pair with SECOND)
SECOND = "#2f6fc4"       # second definition of a pair
HOSP = "#4b57c9"         # hospital or death records (validated pair with PRIM)
PRIM = "#12998c"         # primary care

MONO = "Menlo, 'SF Mono', Consolas, 'Liberation Mono', monospace"
SERIF = "'Iowan Old Style', 'Palatino Linotype', Palatino, Georgia, serif"
MIN_PT = 7.0


def esc(s):
    return html.escape(str(s), quote=True)


def fmt(v, d=2, sign=False):
    """Numbers with an ASCII hyphen for minus, and no '-0.00'."""
    s = f"{v:+.{d}f}" if sign else f"{v:.{d}f}"
    if float(s) == 0:
        s = f"{0:.{d}f}"
    return s


class Figure:
    def __init__(self, width, height):
        self.w, self.h = float(width), float(height)
        self.parts = [f'<rect width="{self.w}" height="{self.h}" fill="{GROUND}"/>']
        self.marks = []          # bounding geometry of data marks: ("seg", x1, y1, x2, y2, halfwidth) or ("box", x0, y0, x1, y1)
        self.texts = []          # (string, x, y, allowed_on_marks)
        self.keyed = False
        self.no_key_reason = None
        self.layout_problems = []   # found while drawing, e.g. a sentence wider than its panel
        self.chrome = []            # part ranges (headline, reading panel, source) the paper version leaves to the caption
        self.text_boxes = []        # (part index, top, bottom) of every text, for cropping the paper version
        self.header_bottom = 0.0
        self.reading_top = None

    # ------------------------------------------------------------------ primitives
    def text(self, x, y, s, size=7.0, fill=INK_2, anchor="start", family=None, weight=None, spacing=None,
             italic=False, on_mark=False, opacity=1.0, max_w=None, rotate=None):
        a = [f'x="{x:.2f}"', f'y="{y:.2f}"', f'font-size="{size}"', f'fill="{fill}"', f'text-anchor="{anchor}"',
             'xml:space="preserve"']
        if family:
            a.append(f'font-family="{family}"')
        if weight:
            a.append(f'font-weight="{weight}"')
        if spacing:
            a.append(f'letter-spacing="{spacing}"')
        if italic:
            a.append('font-style="italic"')
        if opacity != 1.0:
            a.append(f'opacity="{opacity:.2f}"')
        if rotate:
            a.append(f'transform="rotate({rotate} {x:.2f} {y:.2f})"')
        self.parts.append(f'<text {" ".join(a)}>{esc(s)}</text>')
        self.texts.append((str(s), x, y, on_mark))
        part_index = len(self.parts) - 1
        # text a PDF viewer clips at the page edge never reaches the extracted text, so bound it here as well
        width = text_width(s, size, serif=(family == SERIF)) + (spacing or 0) * len(str(s))
        x0 = x if anchor == "start" else (x - width / 2 if anchor == "middle" else x - width)
        if rotate == -90:          # reads upward: the text runs from y up to y - width
            bx0, by0, bx1, by1 = x - 0.75 * size, y - width, x + 0.25 * size, y
        else:
            bx0, by0, bx1, by1 = x0, y - 0.75 * size, x0 + width, y + 0.25 * size
        self.text_boxes.append((part_index, by0, by1))
        if bx0 < -0.5 or bx1 > self.w + 0.5 or by0 < -0.5 or by1 > self.h + 0.5:
            self.layout_problems.append(f"text off the canvas: {str(s)[:40]!r}")
        if max_w is not None and width > max_w + 0.5:
            self.layout_problems.append(f"text wider than its space by {width - max_w:.0f} pt: {str(s)[:40]!r}")

    def line(self, x1, y1, x2, y2, stroke=GRID, width=0.6, dash=None, cap="butt", mark=False, opacity=1.0):
        d = f' stroke-dasharray="{dash}"' if dash else ""
        self.parts.append(f'<line x1="{x1:.2f}" y1="{y1:.2f}" x2="{x2:.2f}" y2="{y2:.2f}" stroke="{stroke}" '
                          f'stroke-width="{width}" stroke-linecap="{cap}" opacity="{opacity:.2f}"{d}/>')
        if mark:
            self.marks.append(("seg", x1, y1, x2, y2, width / 2))

    def rect(self, x, y, w, h, fill, rx=0, stroke=None, sw=0.6, mark=False, opacity=1.0):
        s = f' stroke="{stroke}" stroke-width="{sw}"' if stroke else ""
        self.parts.append(f'<rect x="{x:.2f}" y="{y:.2f}" width="{max(w, 0):.2f}" height="{max(h, 0):.2f}" '
                          f'rx="{rx}" fill="{fill}" opacity="{opacity:.2f}"{s}/>')
        if mark:
            self.marks.append(("box", x, y, x + w, y + h))

    def circle(self, cx, cy, r, fill, stroke=None, sw=0.8, mark=True, opacity=1.0):
        s = f' stroke="{stroke}" stroke-width="{sw}"' if stroke else ""
        self.parts.append(f'<circle cx="{cx:.2f}" cy="{cy:.2f}" r="{max(r, 0):.2f}" fill="{fill}" '
                          f'opacity="{opacity:.2f}"{s}/>')
        if mark:
            self.marks.append(("box", cx - r, cy - r, cx + r, cy + r))

    def path(self, d, stroke=INK, width=1.0, fill="none", dash=None, cap="round", join="round", mark_points=None,
             opacity=1.0):
        da = f' stroke-dasharray="{dash}"' if dash else ""
        self.parts.append(f'<path d="{d}" fill="{fill}" stroke="{stroke}" stroke-width="{width}" '
                          f'stroke-linecap="{cap}" stroke-linejoin="{join}" opacity="{opacity:.2f}"{da}/>')
        if mark_points:
            for (x1, y1), (x2, y2) in zip(mark_points[:-1], mark_points[1:]):
                self.marks.append(("seg", x1, y1, x2, y2, width / 2))

    def polyline(self, pts, stroke=INK, width=1.0, dash=None, mark=True, opacity=1.0):
        d = "M " + " L ".join(f"{x:.2f},{y:.2f}" for x, y in pts)
        self.path(d, stroke=stroke, width=width, dash=dash, mark_points=pts if mark else None, opacity=opacity)

    def polygon(self, pts, fill, stroke="none", width=0.6, mark=False, opacity=1.0):
        d = "M " + " L ".join(f"{x:.2f},{y:.2f}" for x, y in pts) + " Z"
        self.path(d, stroke=stroke, width=width, fill=fill, opacity=opacity)
        if mark:
            xs, ys = [p[0] for p in pts], [p[1] for p in pts]
            self.marks.append(("box", min(xs), min(ys), max(xs), max(ys)))

    # ------------------------------------------------------------------ house furniture
    def header(self, kicker, title, subtitle=None, x=14.0, size=None):
        i0 = len(self.parts)
        size = size or (13.0 if self.w > 300 else 11.0)
        self.text(x, 17, kicker.upper(), 7.0, DIM, spacing=1.4)
        y = 17 + size + 5
        for i, t in enumerate(title if isinstance(title, (list, tuple)) else [title]):
            self.text(x, y + i * (size + 2), t, size, INK, family=SERIF)
        y += (len(title) - 1) * (size + 2) if isinstance(title, (list, tuple)) else 0
        if subtitle:
            for i, t in enumerate(subtitle if isinstance(subtitle, (list, tuple)) else [subtitle]):
                y += 10.5
                self.text(x, y, t, 7.2, DIM)
        self.chrome.append((i0, len(self.parts)))
        self.header_bottom = y + 6
        return y + 6

    def reading(self, x, y, w, lines, title="READING THIS", strong=()):
        i0 = len(self.parts)
        self.reading_top = y
        h = 14 + 10 * len(lines) + 4
        self.rect(x, y, w, h, PANEL, rx=2)
        self.text(x + 8, y + 11, title, 7.0, DIM, spacing=1.2)
        for i, s in enumerate(lines):
            if text_width(s, 7.2) > w - 16:
                self.layout_problems.append(f"reading line wider than its panel by {text_width(s, 7.2) - w + 16:.0f} pt: {s[:40]!r}")
            self.text(x + 8, y + 23 + i * 10, s, 7.2, INK if i in strong else INK_2)
        self.chrome.append((i0, len(self.parts)))
        return y + h

    def key_row(self, x, y, entries, gap=10.0, size=7.0):
        """entries: (kind, colour, label); kind is dot, hollow, square, line, dash, band. Spaced by label width."""
        self.keyed = True
        for kind, colour, label in entries:
            if kind == "dot":
                self.circle(x + 3, y - 2.5, 2.8, colour, mark=False)
            elif kind == "hollow":
                self.circle(x + 3, y - 2.5, 2.5, GROUND, stroke=colour, sw=1.0, mark=False)
            elif kind == "square":
                self.rect(x, y - 6, 6, 6, colour)
            elif kind == "line":
                self.line(x - 1, y - 2.5, x + 8, y - 2.5, colour, 1.6, cap="round")
            elif kind == "dash":
                self.line(x - 1, y - 2.5, x + 8, y - 2.5, colour, 1.2, dash="2 1.6")
            elif kind == "band":
                self.rect(x - 1, y - 5.5, 9, 6, colour, rx=1)
            lx = x + 11
            self.text(lx, y, label, size, INK_2)
            x = lx + text_width(label, size) + gap
        return x

    def key_traits(self, x, y, traits=FIVE, kind="dot"):
        return self.key_row(x, y, [(kind, TRAIT[t], TRAIT_NAME[t]) for t in traits])

    def source(self, s):
        i0 = len(self.parts)
        self.text(self.w - 10, self.h - 6, s, 7.0, DIM, anchor="end")
        self.chrome.append((i0, len(self.parts)))

    def no_key(self, reason):
        self.no_key_reason = reason

    # ------------------------------------------------------------------ output
    def svg(self):
        head = (f'<svg xmlns="http://www.w3.org/2000/svg" width="{self.w}pt" height="{self.h}pt" '
                f'viewBox="0 0 {self.w} {self.h}" font-family="{MONO}">')
        return head + "\n" + "\n".join(self.parts) + "\n</svg>\n"

    def paper_svg(self):
        """The same figure without its headline, reading panel and source line, cropped to what is left.

        In the paper the caption states the finding and the source, so the figure keeps only its key, labels
        and data. Returns the svg, the crop offset and the new height."""
        drop = {i for a, b in self.chrome for i in range(a, b)}
        boxes = [(t_, b_) for i, t_, b_ in self.text_boxes if i not in drop]
        top = max(0.0, min([self.header_bottom] + [t_ for t_, _ in boxes]) - 4)
        bottom = (self.reading_top - 4) if self.reading_top is not None else max(b_ for _, b_ in boxes) + 6
        bottom = max(bottom, max(b_ for _, b_ in boxes) + 4)
        h = bottom - top
        body = [p for i, p in enumerate(self.parts) if i and i not in drop]
        svg = (f'<svg xmlns="http://www.w3.org/2000/svg" width="{self.w}pt" height="{h:.2f}pt" '
               f'viewBox="0 0 {self.w} {h:.2f}" font-family="{MONO}">\n'
               f'<rect width="{self.w}" height="{h:.2f}" fill="{GROUND}"/>\n<g transform="translate(0,{-top:.2f})">\n'
               + "\n".join(body) + "\n</g>\n</svg>\n")
        return svg, top, h

    def save(self, name, check=True, paper=True, png=True):
        """paper=False for a figure drawn only for the paper (a single-column layout); it is saved as given."""
        svg = self.svg()
        (FIG_DIR / f"{name}.svg").write_text(svg)
        cairosvg.svg2pdf(bytestring=svg.encode(), write_to=str(FIG_DIR / f"{name}.pdf"))
        if png:
            cairosvg.svg2png(bytestring=svg.encode(), write_to=str(FIG_DIR / f"{name}.png"),
                             output_width=int(round(self.w * 600 / 72)))
        problems = self.check(FIG_DIR / f"{name}.pdf") if check else []
        if not paper:
            print(f"[{name}] {'CHECKS OK' if not problems else f'PROBLEMS ({len(problems)})'}  {self.w:.0f} x {self.h:.0f} pt")
            for p in problems:
                print("   ", p)
            return problems
        psvg, top, ph = self.paper_svg()
        (FIG_DIR / f"{name}_paper.svg").write_text(psvg)
        cairosvg.svg2pdf(bytestring=psvg.encode(), write_to=str(FIG_DIR / f"{name}_paper.pdf"))
        paper_problems = self.check(FIG_DIR / f"{name}_paper.pdf", dy=top, h=ph, layout=False) if check else []
        problems += [f"paper version: {p}" for p in paper_problems]
        print(f"[{name}] {'CHECKS OK' if not problems else f'PROBLEMS ({len(problems)})'}  {self.w:.0f} x {self.h:.0f} pt, paper version {ph:.0f} pt tall")
        for p in problems:
            print("   ", p)
        return problems

    def check(self, pdf_path, dy=0.0, h=None, layout=True):
        """dy: how far the rendered page was cropped from the top (the paper version); h: that page's height."""
        h = self.h if h is None else h
        problems = list(self.layout_problems) if layout else []
        if layout and not self.keyed and not self.no_key_reason:
            problems.append("no key: declare one with key_row/key_traits, or no_key('why')")
        page = pymupdf.open(pdf_path)[0]
        lines = []
        for b in page.get_text("dict")["blocks"]:
            for ln in b.get("lines", []):
                spans = [s for s in ln["spans"] if s["text"].strip()]
                if not spans:
                    continue
                txt = "".join(s["text"] for s in ln["spans"]).strip()
                x0 = min(s["bbox"][0] for s in spans); x1 = max(s["bbox"][2] for s in spans)
                # glyph box, not the line box: trim the font's ascender and descender padding
                size = max(s["size"] for s in spans)
                if abs(ln["dir"][0]) < 0.5:     # vertical text: the span box is already the glyph column
                    y0_ = min(s["bbox"][1] for s in spans); y1_ = max(s["bbox"][3] for s in spans)
                    ox = spans[0]["origin"][0]
                    lines.append((txt, (ox - 0.72 * size, y0_, ox + 0.18 * size, y1_), min(s["size"] for s in spans)))
                    continue
                base = spans[0]["origin"][1]
                lines.append((txt, (x0, base - 0.72 * size, x1, base + 0.18 * size), min(s["size"] for s in spans)))
        for txt, (x0, y0, x1, y1), size in lines:
            if size < MIN_PT - 0.05:
                problems.append(f"text below {MIN_PT:g} pt ({size:.1f}): {txt[:40]!r}")
            if x0 < -0.5 or y0 < -0.5 or x1 > self.w + 0.5 or y1 > h + 0.5:
                problems.append(f"text off the canvas: {txt[:40]!r}")
        lines = [(txt, (x0, y0 + dy, x1, y1 + dy), size) for txt, (x0, y0, x1, y1), size in lines]   # back to drawing coordinates
        for i in range(len(lines)):
            a = lines[i][1]
            for j in range(i + 1, len(lines)):
                b = lines[j][1]
                if a[0] < b[2] - 0.3 and b[0] < a[2] - 0.3 and a[1] < b[3] - 0.3 and b[1] < a[3] - 0.3:
                    problems.append(f"text overlaps text: {lines[i][0][:30]!r} / {lines[j][0][:30]!r}")
        norm = lambda s: "".join(s.split())
        allowed = [(norm(s), x, y) for s, x, y, ok in self.texts if ok]
        for txt, box, _ in lines:
            key = norm(txt)
            if any(key and (key in s or s in key) and abs(box[3] - 0.18 * 7 - y) < 4 for s, x, y in allowed):
                continue
            for m in self.marks:
                if _hits(m, box):
                    problems.append(f"text on a data mark: {txt[:40]!r}")
                    break
        return problems


def _hits(mark, box, shrink=0.4):
    x0, y0, x1, y1 = box[0] + shrink, box[1] + shrink, box[2] - shrink, box[3] - shrink
    if mark[0] == "box":
        _, a0, b0, a1, b1 = mark
        return a0 < x1 and x0 < a1 and b0 < y1 and y0 < b1
    _, ax, ay, bx, by, hw = mark
    x0, y0, x1, y1 = x0 - hw, y0 - hw, x1 + hw, y1 + hw
    # Liang-Barsky clip of the segment against the text box
    dx, dy = bx - ax, by - ay
    t0, t1 = 0.0, 1.0
    for p, q in ((-dx, ax - x0), (dx, x1 - ax), (-dy, ay - y0), (dy, y1 - ay)):
        if p == 0:
            if q < 0:
                return False
        else:
            t = q / p
            if p < 0:
                t0 = max(t0, t)
            else:
                t1 = min(t1, t)
            if t0 > t1:
                return False
    return True


_WIDTH_CACHE = {}


def text_width(s, size=7.0, serif=False):
    """Width of a monospace label in points (Menlo advance is 0.602 em); serif widths are measured by rendering."""
    if not serif:
        return 0.602 * size * len(str(s))
    key = (s, size)
    if key not in _WIDTH_CACHE:
        svg = (f'<svg xmlns="http://www.w3.org/2000/svg" width="2000pt" height="40pt" viewBox="0 0 2000 40">'
               f'<text x="0" y="30" font-size="{size}" font-family="{SERIF}">{esc(s)}</text></svg>')
        doc = pymupdf.open("pdf", cairosvg.svg2pdf(bytestring=svg.encode()))
        words = doc[0].get_text("words")
        _WIDTH_CACHE[key] = max(w[2] for w in words) if words else 0.0
    return _WIDTH_CACHE[key]


def scale(lo, hi, a, b):
    """Linear map from data [lo, hi] to page [a, b]."""
    return lambda v: a + (v - lo) / (hi - lo) * (b - a)
