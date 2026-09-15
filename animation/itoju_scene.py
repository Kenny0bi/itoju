"""itoju: LD score regression, the one idea the study rests on.

Render:  PATH=~/Library/TinyTeX/bin/universal-darwin:$PATH manim -qh animation/itoju_scene.py LDScoreRegression

About 45 seconds. Every number on screen comes from animation/ldsc_numbers.json, written by
animation/export_ldsc_numbers.py from this study's summary statistics and LDSC output.

  1. A SNP's LD score counts how much of its neighbourhood it tags. For two GWAS, the product of a
     SNP's z-scores grows with its LD score:  E[z1 z2 | l] = (sqrt(N1 N2) rho_g / M) l + rho N_s / sqrt(N1 N2).
  2. Autism x constipation (FinnGen): 950,551 SNPs in 20 LD score bins. The slope is the genetic
     covariance; the intercept sits at zero because no one is in both studies.
  3. EHR depression x constipation: a far steeper slope, and the intercept lifts off zero, because
     FinnGen participants are inside the depression GWAS. Shared people move the intercept, not the slope.
  4. rg = genetic covariance / sqrt(h2_1 h2_2): 0.11 for autism, 0.52 for EHR depression.
"""
import json
from pathlib import Path

import numpy as np
from manim import (
    DOWN, LEFT, RIGHT, UP, Axes, Create, DashedLine, Dot, FadeIn, FadeOut, LaggedStart, Line,
    MathTex, Scene, Text, VGroup, Write, config, Indicate,
)

NUM = json.loads((Path(__file__).resolve().parent / "ldsc_numbers.json").read_text())
P = {p["trait"]: p for p in NUM["pairs"]}

VIOLET = "#6a48c4"     # autism, as in every figure
TURMERIC = "#b98500"   # depression; darker step than print so it reads on white video
INK = "#1b1b1b"
MUTED = "#8a8882"
config.background_color = "#ffffff"
Text.set_default(font="Helvetica Neue")


def fmt(v, d=3):
    # LaTeX typesets "-" as a true minus and rejects U+2212, so numbers keep the ASCII hyphen.
    return f"{v:.{d}f}"


class LDScoreRegression(Scene):
    def panel(self, pair, color, x_center, y_max, y_ticks, title):
        ax = Axes(x_range=[0, 90, 30], y_range=[-y_max * 0.2, y_max, y_ticks[1] - y_ticks[0]], x_length=5.2, y_length=3.6,
                  axis_config={"color": MUTED, "stroke_width": 2, "include_ticks": True}, tips=False).move_to([x_center, -0.55, 0])
        # tick labels in the same sans face as everything else (Axes' own labels are LaTeX serif)
        ticks = VGroup(*[Text(f"{v:g}", font_size=16, color=MUTED).next_to(ax.c2p(v, -y_max * 0.2), DOWN, buff=0.12)
                         for v in (0, 30, 60, 90)],
                       *[Text(f"{v:g}", font_size=16, color=MUTED).next_to(ax.c2p(0, v), LEFT, buff=0.12) for v in y_ticks])
        # the axis line Manim draws at y = 0 is hidden: the dashed zero line marks zero, and a baseline
        # at the bottom carries the tick labels, so "LD score" never collides with them
        ax.x_axis.set_opacity(0)
        base = Line(ax.c2p(0, -y_max * 0.2), ax.c2p(90, -y_max * 0.2), color=MUTED, stroke_width=2)
        ticks.add(base)
        xl = Text("LD score", font_size=20, color=INK).next_to(ax.c2p(45, -y_max * 0.2), DOWN, buff=0.5)
        yl = MathTex(r"\overline{z_1 z_2}", font_size=30, color=INK).next_to(ax.y_axis, UP, buff=0.12)
        head = Text(title, font_size=22, color=INK).next_to(ax, UP, buff=0.55)
        dots = VGroup(*[Dot(ax.c2p(l, z), radius=0.055, color=color) for l, z in zip(pair["bins_l2"], pair["bins_zz"])])
        line = Line(ax.c2p(0, pair["intercept"]), ax.c2p(90, pair["intercept"] + pair["slope"] * 90),
                    color=color, stroke_width=4)
        zero = DashedLine(ax.c2p(0, 0), ax.c2p(90, 0), color=MUTED, stroke_width=1.5, dash_length=0.08)
        return VGroup(ax, ticks), ax, xl, yl, head, dots, line, zero

    def construct(self):
        title = Text("One line tells two stories", font_size=40, color=INK).move_to(UP * 1.6)
        self.play(Write(title), run_time=1.0)
        eq = MathTex(r"\mathbb{E}[\,z_1 z_2 \mid \ell\,]", r"=", r"\frac{\sqrt{N_1 N_2}\,\rho_g}{M}", r"\,\ell",
                     r"+", r"\frac{\rho\, N_s}{\sqrt{N_1 N_2}}", font_size=48, color=INK).next_to(title, DOWN, buff=0.5)
        eq[2].set_color(VIOLET)
        eq[5].set_color(TURMERIC)
        self.play(Write(eq), run_time=2.0)
        k1 = Text("slope: genetic covariance", font_size=24, color=INK)
        k2 = Text("intercept: people in both studies", font_size=24, color=INK)
        keys = VGroup(k1, k2).arrange(RIGHT, buff=1.0).next_to(eq, DOWN, buff=0.6)
        u1 = Line(k1.get_corner(DOWN + LEFT), k1.get_corner(DOWN + RIGHT), color=VIOLET, stroke_width=3).shift(DOWN * 0.1)
        u2 = Line(k2.get_corner(DOWN + LEFT), k2.get_corner(DOWN + RIGHT), color=TURMERIC, stroke_width=3).shift(DOWN * 0.1)
        ell = Text("l = LD score: how much of its neighbourhood a SNP tags", font_size=20, color=MUTED).next_to(keys, DOWN, buff=0.5)
        self.play(FadeIn(keys), Create(u1), Create(u2), FadeIn(ell), run_time=1.2)
        self.wait(2.2)
        self.play(FadeOut(VGroup(keys, u1, u2, ell, title)), eq.animate.scale(0.62).to_edge(UP, buff=0.3), run_time=0.9)

        # 2. autism x constipation
        a = P["ASD"]
        grp1, ax1, xl1, yl1, h1, d1, l1, z1 = self.panel(a, VIOLET, -3.4, 0.08, (0, 0.04, 0.08), "Autism × constipation (FinnGen)")
        self.play(Create(grp1), FadeIn(xl1), FadeIn(yl1), FadeIn(h1), run_time=1.2)
        # font 22: at 18 pt Pango drops some word spaces in this face ("ofabout")
        note = Text(f"Each dot averages about {round(a['n_snps'] / 20, -2):,.0f} SNPs with similar LD scores, {a['n_snps']:,} in all",
                    font_size=22, color=MUTED).next_to(eq, DOWN, buff=0.15)
        self.play(FadeIn(note), LaggedStart(*[FadeIn(d, scale=0.5) for d in d1], lag_ratio=0.08), Create(z1), run_time=2.2)
        self.play(Create(l1), run_time=1.2)
        s1 = Text(f"slope → covariance {a['gencov_obs']:.4f}", font_size=20, color=VIOLET)
        i1 = Text(f"intercept {a['intercept']:.3f}: no one in both studies", font_size=20, color=INK)
        g1 = VGroup(s1, i1).arrange(DOWN, aligned_edge=LEFT, buff=0.12).move_to(ax1.c2p(48, 0.074))
        self.play(FadeIn(g1), Indicate(z1, color=INK), run_time=1.3)
        self.wait(1.4)

        # 3. EHR depression x constipation
        m = P["MDD_EHR"]
        # y axis runs to 1.2 so the steep line leaves the upper left free for the text block
        grp2, ax2, xl2, yl2, h2, d2, l2, z2 = self.panel(m, TURMERIC, 3.4, 1.2, (0, 0.4, 0.8, 1.2), "Depression (EHR) × constipation")
        self.play(Create(grp2), FadeIn(xl2), FadeIn(yl2), FadeIn(h2), run_time=1.2)
        self.play(LaggedStart(*[FadeIn(d, scale=0.5) for d in d2], lag_ratio=0.08), Create(z2), run_time=2.0)
        self.play(Create(l2), run_time=1.2)
        s2 = Text(f"slope → covariance {m['gencov_obs']:.4f}", font_size=20, color=TURMERIC)
        i2 = Text(f"intercept {m['intercept']:.3f}:", font_size=20, color=INK)
        w2 = Text("FinnGen people in both studies", font_size=20, color=INK)
        g2 = VGroup(s2, i2, w2).arrange(DOWN, aligned_edge=LEFT, buff=0.12).next_to(ax2.c2p(5, 1.17), DOWN + RIGHT, buff=0)
        lift = Line(ax2.c2p(0, 0), ax2.c2p(0, m["intercept"]), color=INK, stroke_width=9)
        lift_dot = Dot(ax2.c2p(0, m["intercept"]), radius=0.09, color=INK)
        self.play(Create(lift), FadeIn(lift_dot, scale=0.5), FadeIn(g2), run_time=1.4)
        self.wait(1.0)
        shared = Text("Shared people lift the intercept. The slope keeps the genetics.", font_size=24, color=INK)
        shared.to_edge(DOWN, buff=0.25)
        self.play(FadeIn(shared), run_time=1.0)
        self.wait(2.2)

        # 4. from covariance to correlation
        everything = VGroup(grp1, xl1, yl1, h1, d1, l1, z1, note, g1, grp2, xl2, yl2, h2, d2, l2, z2, g2, lift, lift_dot, shared)
        self.play(FadeOut(everything), run_time=0.9)
        rg_eq = MathTex(r"r_g", r"=", r"\frac{\rho_g}{\sqrt{h^2_1\, h^2_2}}", font_size=56, color=INK).move_to(UP * 1.0)
        self.play(Write(rg_eq), run_time=1.4)
        r1 = Text(f"autism × constipation:  rg = {a['rg']:.2f} ± {a['rg_se']:.2f}", font_size=30, color=VIOLET)
        r2 = Text(f"EHR depression × constipation:  rg = {m['rg']:.2f} ± {m['rg_se']:.2f}", font_size=30, color=TURMERIC)
        rows = VGroup(r1, r2).arrange(DOWN, buff=0.35).next_to(rg_eq, DOWN, buff=0.7)
        self.play(LaggedStart(FadeIn(r1), FadeIn(r2), lag_ratio=0.6), run_time=2.0)
        close = Text("itoju asks how this number moves when only the definition changes.",
                     font_size=24, color=INK).to_edge(DOWN, buff=0.7)
        self.play(FadeIn(close), run_time=1.0)
        self.wait(3.0)
