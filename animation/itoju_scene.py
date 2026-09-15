"""itoju: why a small shift in genetic correlation between two definitions can still be certain.

The one idea the study's test rests on, drawn with its real numbers (animation/block_numbers.json,
written by animation/export_block_numbers.py from the LDSC runs the paper uses):

  1. Depression's genetic correlation with sleep apnoea is 0.353 under hospital records and 0.420
     under any sleep disorder. On their own, the two 95% intervals overlap.
  2. LDSC measures its noise by leaving out one of 200 genome blocks at a time. On a common SNP set,
     block k is the same stretch of genome under both definitions, so every block gives a pair of
     pseudovalues, and the pairs lie along the diagonal (r = 0.96).
  3. Subtract block by block and the shared noise cancels: the standard error of the difference is
     0.007. Pair each block with a random other block, as if the definitions shared no people, and
     it is 0.035.
  4. The shift of +0.068 set against its own noise.

Render (TinyTeX on PATH for the one formula):
  PATH=~/Library/TinyTeX/bin/universal-darwin:$PATH manim -qh animation/itoju_scene.py DefinitionShift
"""
import json
from pathlib import Path

import numpy as np
from manim import (DOWN, LEFT, ORIGIN, RIGHT, UP, Create, DashedLine, Dot, FadeIn, FadeOut, GrowFromCenter,
                   LaggedStart, Line, MathTex, Rectangle, Scene, Text, VGroup, config, linear)

B = json.loads((Path(__file__).resolve().parent / "block_numbers.json").read_text())

GROUND = "#EEF3F1"
PANEL = "#E2EAE7"
INK = "#14232B"
INK_2 = "#3A4A51"
DIM = "#5E6E72"
RULE = "#A9B8B2"
NOISE = "#CDD8D4"
MDD = "#e8a300"
FIRST = "#c9761a"
SECOND = "#2f6fc4"
SERIF = "Charter"
MONO = "Menlo"
config.background_color = GROUND


def T(s, size=24, color=INK, font=MONO):
    """Text at its final size (drawing it larger and scaling down makes Pango wrap long lines)."""
    return Text(s, font=font, font_size=size, color=color)


def inside(m, margin=0.3):
    """Shift a mobject sideways until it sits inside the frame, so no label is ever cut at an edge."""
    half = config.frame_width / 2 - margin
    if m.get_right()[0] > half:
        m.shift(LEFT * (m.get_right()[0] - half))
    if m.get_left()[0] < -half:
        m.shift(RIGHT * (-half - m.get_left()[0]))
    return m


def swarm(values, x_of, y0, r):
    """Beeswarm positions: each dot as close to the baseline as it can sit without touching another."""
    d = 2 * r * 1.02
    placed, pos = [], [None] * len(values)
    for i in np.argsort(np.abs(values - np.median(values))):
        x = x_of(values[i])
        k = 0
        while True:
            off = ((k + 1) // 2) * d * (1 if k % 2 else -1)
            if all((x - px) ** 2 + (off - py) ** 2 >= d * d for px, py in placed):
                break
            k += 1
        placed.append((x, off))
        pos[i] = np.array([x, y0 + off, 0.0])
    return pos


class DefinitionShift(Scene):
    def construct(self):
        rg1, rg2, se1, se2 = B["rg1"], B["rg2"], B["se1"], B["se2"]
        shift, se_pair, se_broken = B["shift"], B["se_paired"], B["se_broken"]
        dev1, dev2 = np.array(B["dev1"]), np.array(B["dev2"])
        n = B["n_blocks"]

        # ------------------------------------------------------------------ title
        kicker = T("ITOJU", 22, DIM)
        title = VGroup(T("Two definitions of one condition,", 46, INK, SERIF),
                       T("measured on the same genome", 46, INK, SERIF)).arrange(DOWN, buff=0.18)
        sub = T("why a small shift in genetic correlation can still be certain", 22, DIM)
        VGroup(kicker, title, sub).arrange(DOWN, buff=0.4)
        self.play(FadeIn(kicker), FadeIn(title, shift=UP * 0.2), run_time=1.2)
        self.play(FadeIn(sub), run_time=0.7)
        self.wait(1.6)
        self.play(FadeOut(VGroup(kicker, title, sub)), run_time=0.6)

        # ------------------------------------------------------------------ 1. two estimates
        head = T("Depression and sleep apnoea: one genetic correlation, two definitions", 28, INK, SERIF)
        head.to_edge(UP, buff=0.5)
        lo_rg, hi_rg = 0.26, 0.50
        X = lambda v: -5.2 + (v - lo_rg) / (hi_rg - lo_rg) * 10.4
        axis_y = -1.3
        axis = Line([X(lo_rg), axis_y, 0], [X(hi_rg), axis_y, 0], color=RULE, stroke_width=2)
        ticks = VGroup()
        for v in (0.30, 0.35, 0.40, 0.45):
            ticks.add(Line([X(v), axis_y, 0], [X(v), axis_y - 0.1, 0], color=RULE, stroke_width=2))
            ticks.add(T(f"{v:.2f}", 21, DIM).move_to([X(v), axis_y - 0.35, 0]))
        axis_lab = T("genetic correlation with depression", 21, DIM).move_to([0, axis_y - 0.8, 0])
        self.play(FadeIn(head), Create(axis), FadeIn(ticks), FadeIn(axis_lab), run_time=1.2)

        rows = VGroup()
        for y, rg, se, col, name in ((1.0, rg1, se1, FIRST, "sleep apnoea, hospital records"),
                                     (-0.1, rg2, se2, SECOND, "any sleep disorder")):
            band = Line([X(rg - 1.96 * se), y, 0], [X(rg + 1.96 * se), y, 0], color=col, stroke_width=14,
                        stroke_opacity=0.28)
            dot = Dot([X(rg), y, 0], radius=0.11, color=col)
            lab = T(name, 22, INK_2).move_to([X(rg - 1.96 * se), y + 0.42, 0], aligned_edge=LEFT)
            val = T(f"{rg:.3f}", 22, INK).next_to([X(rg + 1.96 * se), y, 0], RIGHT, buff=0.3)
            rows.add(VGroup(band, dot, lab, val))
            self.play(GrowFromCenter(band), FadeIn(dot, scale=0.5), FadeIn(lab), FadeIn(val), run_time=1.0)
        note = T("On their own, the two 95% intervals overlap. Is the shift real?", 24, INK).to_edge(DOWN, buff=0.45)
        self.play(FadeIn(note), run_time=0.7)
        self.wait(1.8)
        self.play(FadeOut(VGroup(rows, axis, ticks, axis_lab, note, head)), run_time=0.6)

        # ------------------------------------------------------------------ 2. 200 blocks, paired
        head2 = T("LDSC measures its noise by leaving out one genome block at a time", 28, INK, SERIF)
        head2.to_edge(UP, buff=0.5)
        strip_y, strip_w = 2.45, 11.6
        cw = strip_w / n
        chrom = B.get("block_chrom", [1] * n)
        strip = VGroup(*[Rectangle(width=cw * 0.8, height=0.34, fill_color="#C9D5D0" if chrom[k] % 2 else "#E0E8E4",
                                   fill_opacity=1, stroke_width=0).move_to([-strip_w / 2 + (k + 0.5) * cw, strip_y, 0])
                         for k in range(n)])
        strip_lab = VGroup(T(f"the genome in {n} blocks, chromosome 1 to 22", 21, DIM)
                           .next_to(strip, DOWN, buff=0.12).align_to(strip, LEFT))
        for c in (1, 6, 12, 22):
            k = chrom.index(c)
            strip_lab.add(T(str(c), 19, DIM).next_to(strip[k], UP, buff=0.08))
        self.play(FadeIn(head2), LaggedStart(*[FadeIn(b) for b in strip], lag_ratio=0.004), FadeIn(strip_lab),
                  run_time=1.3)

        lim, side = 1.2, 4.6
        cx, cy = -3.3, -1.0
        P = lambda a, b: np.array([cx + np.clip(a, -lim, lim) / lim * side / 2,
                                   cy + np.clip(b, -lim, lim) / lim * side / 2, 0])
        frame = Rectangle(width=side, height=side, stroke_color=RULE, stroke_width=1.5).move_to([cx, cy, 0])
        diag = DashedLine(P(-lim, -lim), P(lim, lim), color=RULE, stroke_width=1.8, dash_length=0.08)
        xl = T("hospital records", 21, DIM).next_to(frame, DOWN, buff=0.14)
        yl = T("any sleep disorder", 21, DIM).rotate(np.pi / 2).next_to(frame, LEFT, buff=0.14)
        dl = T("each dot: one block's pull on both estimates", 21, DIM).next_to(frame, UP, buff=0.14).align_to(frame, LEFT)
        self.play(Create(frame), Create(diag), FadeIn(xl), FadeIn(yl), FadeIn(dl), run_time=0.8)

        dots = VGroup(*[Dot(P(a, b), radius=0.04, color=MDD) for a, b in zip(dev1, dev2)])
        sweep = Rectangle(width=cw * 1.6, height=0.5, stroke_color=INK, stroke_width=2.5).move_to(strip[0])
        self.add(sweep)
        self.play(sweep.animate(rate_func=linear).move_to(strip[-1]),
                  LaggedStart(*[FadeIn(d, scale=0.3) for d in dots], lag_ratio=0.02, rate_func=linear),
                  run_time=3.6)
        self.play(FadeOut(sweep), run_time=0.3)
        r_text = VGroup(T(f"blocks agree, r = {B['error_corr']:.2f}", 30, INK, SERIF),
                        T("the same stretch of genome", 22, INK_2),
                        T("pushes both estimates the same way,", 22, INK_2),
                        T("because the definitions share people", 22, INK_2)).arrange(DOWN, aligned_edge=LEFT, buff=0.16)
        r_text.move_to([-0.5, cy + 0.4, 0], aligned_edge=LEFT)
        assert r_text.get_right()[0] < config.frame_width / 2 - 0.3, r_text.get_right()
        self.play(FadeIn(r_text, shift=LEFT * 0.2), run_time=0.9)
        self.wait(2.0)

        # ------------------------------------------------------------------ 3. subtract block by block
        head3 = T("Subtract block by block, and the shared noise cancels", 28, INK, SERIF).to_edge(UP, buff=0.5)
        self.play(FadeOut(VGroup(strip, strip_lab, frame, diag, xl, yl, dl, r_text)), FadeOut(head2), FadeIn(head3),
                  run_time=0.8)
        lo_d, hi_d = -1.6, 1.6
        D = lambda v: -5.4 + (np.clip(v, lo_d, hi_d) - lo_d) / (hi_d - lo_d) * 10.8
        y_pair, y_broken = 0.9, -2.0
        diff = dev2 - dev1
        rng = np.random.default_rng(7)
        broken = dev2[rng.permutation(n)] - dev1
        pos_pair = swarm(diff, D, y_pair, 0.04)
        pos_broken = swarm(broken, D, y_broken, 0.04)
        half_pair = max(abs(p[1] - y_pair) for p in pos_pair) + 0.12
        base_pair = Line([D(lo_d), y_pair - half_pair, 0], [D(hi_d), y_pair - half_pair, 0], color=RULE, stroke_width=1.5)
        base_broken = Line([D(lo_d), y_broken - 0.55, 0], [D(hi_d), y_broken - 0.55, 0], color=RULE, stroke_width=1.5)
        zero = DashedLine([D(0), y_pair + half_pair, 0], [D(0), y_broken - 0.55, 0], color=RULE, stroke_width=1.5,
                          dash_length=0.08)
        self.play(*[d.animate.move_to(p) for d, p in zip(dots, pos_pair)], Create(base_pair), run_time=2.4)
        formula = MathTex(r"\mathrm{SE}(\Delta)=\sqrt{\mathrm{var}\left(\tilde r_{2,k}-\tilde r_{1,k}\right)/n}",
                          color=INK, font_size=34).next_to(head3, DOWN, buff=0.3)
        lab_pair = VGroup(T("blocks paired as they are", 22, INK_2), T(f"SE {se_pair:.3f}", 30, INK, SERIF)) \
            .arrange(DOWN, aligned_edge=LEFT, buff=0.1).move_to([D(hi_d) - 1.2, y_pair + 0.55, 0])
        inside(lab_pair)
        self.play(FadeIn(formula), FadeIn(lab_pair), Create(zero), run_time=0.9)
        self.wait(1.2)
        ghosts = VGroup(*[Dot(p, radius=0.04, color=RULE) for p in pos_pair])
        self.add(ghosts)
        self.play(*[g.animate.move_to(p) for g, p in zip(ghosts, pos_broken)], Create(base_broken), run_time=2.4)
        lab_broken = VGroup(T("each block paired with a random one,", 22, INK_2),
                            T("as if the definitions shared no people", 22, INK_2),
                            T(f"SE {se_broken:.3f}", 30, INK, SERIF)).arrange(DOWN, aligned_edge=LEFT, buff=0.1)
        lab_broken.move_to([D(hi_d) - 1.8, y_broken + 1.0, 0])
        inside(lab_broken)
        self.play(FadeIn(lab_broken), run_time=0.8)
        self.wait(2.2)

        # ------------------------------------------------------------------ 4. the shift against its noise
        head4 = T("The shift, set against its own noise", 28, INK, SERIF).to_edge(UP, buff=0.5)
        self.play(FadeOut(VGroup(dots, ghosts, base_pair, base_broken, zero, formula, lab_pair, lab_broken)),
                  FadeOut(head3), FadeIn(head4), run_time=0.8)
        lo_s, hi_s = -0.10, 0.10
        S = lambda v: -5.2 + (v - lo_s) / (hi_s - lo_s) * 8.0
        ax_y = -2.1
        ax = Line([S(lo_s), ax_y, 0], [S(hi_s), ax_y, 0], color=RULE, stroke_width=2)
        tk = VGroup()
        for v, s in ((-0.1, "-0.10"), (-0.05, "-0.05"), (0, "0"), (0.05, "+0.05"), (0.1, "+0.10")):
            tk.add(Line([S(v), ax_y, 0], [S(v), ax_y - 0.1, 0], color=RULE, stroke_width=2))
            tk.add(T(s, 18, DIM).move_to([S(v), ax_y - 0.35, 0]))
        tk.add(T("shift in rg, any sleep disorder minus hospital records", 21, DIM)
               .move_to([0, ax_y - 0.8, 0]))
        inside(tk[-1])
        assert tk[-1].width < config.frame_width - 0.6, tk[-1].width
        self.play(Create(ax), FadeIn(tk), run_time=0.8)

        inside_broken = abs(shift) <= 1.96 * se_broken
        rows4 = VGroup()
        for y, hw, name, verdict in (
                (0.9, 1.96 * se_pair, "same people, blocks paired", "far outside its noise"),
                (-0.7, 1.96 * se_broken, "if no people were shared",
                 "inside its noise" if inside_broken else "just outside its noise")):
            sleeve = Rectangle(width=S(hw) - S(-hw), height=0.62, fill_color=NOISE, fill_opacity=1, stroke_width=0)
            sleeve.move_to([(S(hw) + S(-hw)) / 2, y, 0])
            zero_tick = Line([S(0), y - 0.45, 0], [S(0), y + 0.45, 0], color=RULE, stroke_width=2)
            stem = Line([S(0), y, 0], [S(shift), y, 0], color=MDD, stroke_width=9)
            head_dot = Dot([S(shift), y, 0], radius=0.12, color=MDD)
            lab = T(name, 22, INK_2).move_to([S(lo_s), y + 0.58, 0], aligned_edge=LEFT)
            ver = T(verdict, 22, INK).next_to([S(shift), y, 0], RIGHT, buff=0.35)
            rows4.add(VGroup(sleeve, zero_tick, stem, head_dot, lab, ver))
            self.play(GrowFromCenter(sleeve), Create(zero_tick), FadeIn(lab), run_time=0.8)
            self.play(Create(stem), FadeIn(head_dot, scale=0.5), run_time=0.8)
            self.play(FadeIn(ver), run_time=0.5)
        mant, expo = f"{B['p']:.0e}".split("e")
        pval = MathTex(rf"\Delta = +{shift:.3f},\quad p = {mant}\times 10^{{{int(expo)}}}", color=INK, font_size=36)
        pval.next_to(head4, DOWN, buff=0.35)
        self.play(FadeIn(pval), run_time=0.7)
        self.wait(2.4)

        # ------------------------------------------------------------------ end
        self.play(FadeOut(VGroup(head4, ax, tk, rows4, pval)), run_time=0.7)
        end = T("itoju", 64, INK, SERIF)
        end2 = T("how a definition in the health record changes the genetics", 22, DIM).next_to(end, DOWN, buff=0.3)
        self.play(FadeIn(end, shift=UP * 0.2), FadeIn(end2), run_time=1.0)
        self.wait(1.8)
