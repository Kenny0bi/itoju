"""Figure: the anatomy of a definition, what each one actually counts.

Every alternative definition is one row. Left block: one tile per ICD-10 code used anywhere in that
condition's definitions, in code order. An indigo tile counts the code in hospital or death records; a teal
corner means the code also counts from primary care. Right block, aligned for every condition: the rules
codes cannot show, each drawn as its own small glyph, in the manner of Dear Data's hand-built marks, with
one shared key at the top:
  pill        counts a drug purchase class (ATC A06A, laxatives)
  card        counts a drug reimbursement right
  slash       the case must not also meet another definition (focal and not generalized)
  bars        the code must be the person's most frequent diagnosis of that kind
  target      the control rule differs from the set's first definition
A faint dot means the rule does not apply. Two definitions with identical code rows can still differ in the
right block, which is exactly where some of them differ. Tiles show ICD-10 only; ICD-9 and ICD-8 rules are
in the endpoint file.

Sources: results/definitions/code_membership.tsv and endpoint_anatomy.tsv (src/02). Every sentence in the
reading panel is computed.
"""
from itertools import combinations

import pandas as pd

import itoju_svg as sv
from itoju_labels import ENDPOINT_LABEL

ROOT = sv.ROOT
SETS = [
    ("Epilepsy", ["G6_EPLEPSY", "FE", "FE_STRICT", "FE_MODE", "GE", "GE_STRICT", "GE_MODE"]),
    ("Sleep apnoea", ["G6_SLEEPAPNO", "G6_SLEEPAPNO_INCLAVO", "SLEEP"]),
    ("Insomnia", ["F5_INSOMNIA", "KRA_PSY_SLEEP_NONORG_EXMORE"]),
    ("Constipation", ["K11_CONSTIPATION", "K11_OTHFUNC"]),
    ("ADHD", ["F5_ADHD", "KRA_PSY_HYPERKIN_EXMORE"]),
]
RULES = [("drug", "drug"), ("reimb", "reimb."), ("excl", "excl."), ("mode", "mode"), ("ctrl", "control")]
LAB_R, TX0, TP = 172.0, 180.0, 11.0      # label right edge, first tile centre, tile pitch
RX0, RP = 372.0, 30.0                     # first rule glyph centre, rule pitch
ROW = 11.0


def glyph(f, kind, x, y, on=True):
    ink = sv.INK_2
    if not on:
        f.circle(x, y, 1.1, sv.GRID, mark=False)
        return
    if kind == "drug":
        f.rect(x - 5, y - 2.3, 10, 4.6, "#FFFFFF", rx=2.3, stroke=ink, sw=0.8)
        f.path(f"M {x:.2f},{y - 2.3:.2f} L {x - 2.7:.2f},{y - 2.3:.2f} A 2.3 2.3 0 0 0 {x - 2.7:.2f},{y + 2.3:.2f} "
               f"L {x:.2f},{y + 2.3:.2f} Z", stroke=ink, width=0.8, fill=ink)
    elif kind == "reimb":
        f.rect(x - 4.7, y - 3.3, 9.4, 6.6, "#FFFFFF", stroke=ink, sw=0.8)
        f.rect(x - 4.7, y - 2.0, 9.4, 1.5, ink)
    elif kind == "excl":
        f.circle(x, y, 3.6, "#FFFFFF", stroke=ink, sw=0.9, mark=False)
        f.line(x - 2.5, y + 2.5, x + 2.5, y - 2.5, ink, 0.9)
    elif kind == "mode":
        for dx, h, fill in ((-3.0, 3.4, False), (0.0, 7.2, True), (3.0, 4.6, False)):
            f.rect(x + dx - 1.1, y + 3.6 - h, 2.2, h, ink if fill else "#FFFFFF", stroke=ink, sw=0.6)
    elif kind == "ctrl":
        f.circle(x, y, 3.8, "#FFFFFF", stroke=ink, sw=0.8, mark=False)
        f.circle(x, y, 1.5, ink, mark=False)


def tile(f, x, y, hosp, prim):
    if hosp:
        f.rect(x - 4.5, y - 4.5, 9, 9, sv.HOSP, rx=0.8)
        if prim:
            f.polygon([(x + 4.5, y - 4.5), (x + 4.5, y + 0.6), (x - 0.6, y - 4.5)], sv.PRIM, stroke=sv.GROUND, width=0.5)
    else:
        f.rect(x - 4.5, y - 4.5, 9, 9, sv.PANEL, rx=0.8, stroke=sv.GRID, sw=0.4)


def main():
    mem = pd.read_csv(ROOT / "results" / "definitions" / "code_membership.tsv", sep="\t")
    ana = pd.read_csv(ROOT / "results" / "definitions" / "endpoint_anatomy.tsv", sep="\t").set_index("endpoint")

    W = sv.DOUBLE
    left = 14.0
    f = sv.Figure(W, 2000.0)          # drawn tall, trimmed to the content at the end
    f.header("itoju  /  anatomy of a definition",
             "Same condition, different people: what each definition actually counts",
             "Every ICD-10 code a definition counts, and the rules no code list can show")
    f.key_row(left, 70, [])
    kx = left
    for hosp, prim, lab in ((True, False, "counted in hospital or death records"), (True, True, "also counted from primary care"),
                            (False, False, "not counted")):
        tile(f, kx + 4.5, 67.5, hosp, prim)
        f.text(kx + 13, 70, lab, 7.0, sv.INK_2)
        kx += 13 + sv.text_width(lab) + 14
    keys = [("drug", "counts a drug purchase"), ("reimb", "counts a reimbursement right"),
            ("excl", "must not meet another definition")]
    keys2 = [("mode", "must be the most frequent code"), ("ctrl", "control rule differs from the first"),
             (None, "rule does not apply")]
    for yk, row in ((84.0, keys), (97.0, keys2)):
        kx = left
        for kind, lab in row:
            glyph(f, kind, kx + 5, yk - 2.5, kind is not None)
            f.text(kx + 13, yk, lab, 7.0, sv.INK_2)
            kx += 13 + sv.text_width(lab) + 14

    hy = 122.0
    f.text(TX0 - 4.5, hy, "ICD-10 CODES", 7.0, sv.DIM, spacing=0.6)
    f.text(RX0 - 5, hy, "RULES", 7.0, sv.DIM, spacing=0.6)
    for j, (_, lab) in enumerate(RULES):
        f.text(RX0 + j * RP, hy + 11, lab, 7.0, sv.INK_2, anchor="middle")
    f.line(left, hy + 16, W - 10, hy + 16, sv.RULE, 0.5)

    y = hy + 22
    rows_by_set, prim_eps, drug_eps = {}, [], []
    for name, defs in SETS:
        codes = sorted(set(mem[mem.endpoint.isin(defs)].icd10))
        codes = [c for c in codes if not (len(c) == 3 and any(o.startswith(c) and len(o) > 3 for o in codes))] or codes
        f.text(left, y + 26, name, 8.0, sv.INK, family=sv.SERIF)
        for k, c in enumerate(codes):
            f.text(TX0 + k * TP + 2.5, y + 27, f"{c[:3]}.{c[3:]}" if len(c) > 3 else c, 7.0, sv.DIM, rotate=-90)
        ref_ctrl = str(ana.loc[defs[0], "control_exclude"]) + str(ana.loc[defs[0], "control_conditions"])
        yy = y + 36
        rows = {}
        for e in defs:
            a = ana.loc[e]
            m = mem[mem.endpoint == e]
            hd = set(m[m.source.isin(["HD", "COD"])].icd10)
            op = set(m[m.source == "OUTPAT"].icd10)
            pattern = []
            for k, c in enumerate(codes):
                h = c in hd or any(x.startswith(c) for x in hd)
                p = c in op or any(x.startswith(c) for x in op)
                tile(f, TX0 + k * TP, yy, h, p)
                pattern.append((h, p))
            sources = str(a.sources).split()
            flags = {"drug": "KELA_ATC" in sources, "reimb": "KELA_REIMB" in sources,
                     "excl": isinstance(a.case_conditions, str) and a.case_conditions.strip() != "",
                     "mode": bool(a.mode_rule),
                     "ctrl": (str(a.control_exclude) + str(a.control_conditions)) != ref_ctrl}
            for j, (kind, _) in enumerate(RULES):
                glyph(f, kind, RX0 + j * RP, yy, flags[kind])
            f.text(LAB_R, yy + 2.5, ENDPOINT_LABEL[e], 7.0, sv.INK_2, anchor="end", max_w=LAB_R - left)
            rows[e] = (tuple(pattern), tuple(flags.values()))
            if any(p for _, p in pattern):
                prim_eps.append((name, e))
            if flags["drug"]:
                drug_eps.append((name, e))
            yy += ROW
        rows_by_set[name] = rows
        y = yy + 4

    same_codes = rule_only = 0
    for rows in rows_by_set.values():
        for a, b in combinations(rows.values(), 2):
            if a[0] == b[0]:
                same_codes += 1
                rule_only += a[1] != b[1]
    n_pairs = sum(len(r) * (len(r) - 1) // 2 for r in rows_by_set.values())
    print(f"  pairs {n_pairs}, identical codes {same_codes}, of which rules differ {rule_only}; primary care {prim_eps}; drug {drug_eps}")
    assert len(prim_eps) == 1 and len(drug_eps) == 1
    def lab(ne):
        l_ = ENDPOINT_LABEL[ne[1]][0].lower() + ENDPOINT_LABEL[ne[1]][1:]
        return l_ if l_.startswith(ne[0].lower()) else f"{ne[0].lower()} {l_}"
    lines = [
        f"{same_codes} of {n_pairs} definition pairs count exactly the same ICD-10 codes; "
        + ("all of them differ only in their rules." if rule_only == same_codes else f"{rule_only} of those differ in their rules."),
        f"Primary care counts only in {lab(prim_eps[0])}; drug purchases only in {lab(drug_eps[0])}.",
        "Tiles show ICD-10 codes only; the ICD-9 and ICD-8 parts of each definition are in the endpoint file.",
    ]
    end = f.reading(left, y + 4, W - left - 10, lines, strong=(0,))
    f.h = end + 20
    f.parts[0] = f'<rect width="{f.w}" height="{f.h}" fill="{sv.GROUND}"/>'
    f.source("FinnGen R12 endpoint definitions; results/definitions")
    f.save("fig02_anatomy")


if __name__ == "__main__":
    main()
