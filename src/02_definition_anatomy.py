"""Turn FinnGen R12 endpoint definitions into comparable code sets and distances.

For every endpoint in the study, this script:
  1. resolves INCLUDE chains, so an endpoint like SLEEP inherits the rules of the endpoints it unions;
  2. expands each registry's ICD-10 patterns against the WHO ICD-10 code universe
     (the codes present in the phecodeX WHO map, plus their three-character parents);
  3. records which data sources can make someone a case: hospital discharge (HD), cause of death
     (COD), primary care (OUTPAT, Avohilmo), drug purchases (KELA_ATC), drug reimbursement
     (KELA_REIMB), and extra case conditions (CONDITIONS) or control rules;
  4. maps the expanded ICD-10 set to phecodeX, the shared standard vocabulary;
  5. computes pairwise definition distances within each condition family.

FinnGen pattern syntax (FinnGen endpoint file format description):
  "|" separates alternatives, each a regular expression matched from the start of the code;
  "%" marks a "mode" rule (the code must be the most common among the related diagnoses);
  "$!$" means the registry is deliberately not used;
  "!NAME" in CONDITIONS means "and not a case of NAME".

Outputs: results/definitions/endpoint_anatomy.tsv, code_membership.tsv, pairwise_distance.tsv
"""
import re
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
ENDPOINTS = ROOT / "data" / "raw" / "endpoints" / "finngen_R12_endpoint_core_noncore_1.0.xlsx"
MANIFEST = ROOT / "data" / "raw" / "finngen_R12_manifest.tsv"
PHECODE_MAP = ROOT / "data" / "ref" / "phecode" / "phecodeX_ICD_WHO_map_flat.csv"
OUT = ROOT / "results" / "definitions"

FAMILIES = {
    "Epilepsy": ["G6_EPLEPSY", "FE", "FE_STRICT", "FE_MODE", "GE", "GE_STRICT", "GE_MODE", "G6_STATUSEPI"],
    "Sleep": ["SLEEP", "G6_SLEEPAPNO", "G6_SLEEPAPNO_INCLAVO", "F5_INSOMNIA",
              "KRA_PSY_SLEEP_NONORG_EXMORE", "G6_SLEEPDISOTH", "F5_SLEEP_NOS"],
    "Bowel": ["K11_CONSTIPATION", "K11_OTHFUNC", "K11_IBS", "K11_FUNCDYSP", "K11_REFLUX"],
    "Bladder": ["N14_NEUROMUSCDYSBLADD", "N14_OTHBLADD"],
    "ADHD": ["F5_ADHD", "KRA_PSY_HYPERKIN_EXMORE"],
    "Intellectual disability": ["F5_MILDRET", "KRA_PSY_MENTALRET_EXMORE"],
    "Autism": ["KRA_PSY_DEVWIDE_EXMORE", "KRA_PSY_AUTISM_EXMORE"],
}

ICD10_COLS = {"HD": "HD_ICD_10", "COD": "COD_ICD_10", "OUTPAT": "OUTPAT_ICD"}
OTHER_SOURCES = {"KELA_ATC": "KELA_ATC", "KELA_REIMB": "KELA_REIMB", "OPER": "OPER_NOM"}


def present(v):
    return pd.notna(v) and str(v).strip() not in ("", "nan", "$!$")


def icd10_universe():
    m = pd.read_csv(PHECODE_MAP, dtype=str, encoding="latin-1")
    codes = set(m["icd"].str.replace(".", "", regex=False))
    codes |= {c[:3] for c in codes}
    return m, sorted(codes)


def split_alternatives(pattern):
    """Split on "|" outside square brackets. FinnGen writes F84[0|5], where the bar sits inside a
    character class, so a plain split would cut the class in half."""
    parts, depth, cur = [], 0, ""
    for ch in pattern:
        if ch == "[":
            depth += 1
        elif ch == "]":
            depth -= 1
        if ch == "|" and depth == 0:
            parts.append(cur)
            cur = ""
        else:
            cur += ch
    parts.append(cur)
    return parts


def expand(pattern, universe):
    """Codes in the universe matched by a FinnGen pattern string."""
    if not present(pattern):
        return set(), False
    out, mode = set(), False
    for alt in split_alternatives(str(pattern)):
        alt = alt.strip()
        if alt.startswith("%"):
            mode, alt = True, alt[1:]
        alt = alt.replace(".", "")
        rx = re.compile("^" + alt)
        out |= {c for c in universe if rx.match(c)}
    return out, mode


def resolve(name, table, universe, seen=None):
    """Case rules for an endpoint, with INCLUDE chains folded in."""
    seen = seen or set()
    if name in seen:
        return None
    seen.add(name)
    r = table.loc[name]
    rule = {"icd10": {s: set() for s in ICD10_COLS}, "sources": set(), "mode": False,
            "atc": set(), "conditions": [], "includes": []}
    for src, col in ICD10_COLS.items():
        codes, mode = expand(r.get(col), universe)
        if codes:
            rule["icd10"][src] |= codes
            rule["sources"].add(src)
        rule["mode"] |= mode
    for src, col in OTHER_SOURCES.items():
        if present(r.get(col)):
            rule["sources"].add(src)
            if src == "KELA_ATC":
                rule["atc"] |= set(str(r[col]).split("|"))
    if present(r.get("CONDITIONS")):
        rule["conditions"].append(str(r["CONDITIONS"]))
    if present(r.get("INCLUDE")):
        for child in str(r["INCLUDE"]).split("|"):
            child = child.strip()
            rule["includes"].append(child)
            sub = resolve(child, table, universe, seen)
            if sub is None:
                continue
            for s in ICD10_COLS:
                rule["icd10"][s] |= sub["icd10"][s]
            rule["sources"] |= sub["sources"]
            rule["mode"] |= sub["mode"]
            rule["atc"] |= sub["atc"]
            rule["conditions"] += sub["conditions"]
    return rule


def jaccard(a, b):
    return len(a & b) / len(a | b) if (a | b) else np.nan


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    table = pd.read_excel(ENDPOINTS).set_index("NAME")
    man = pd.read_csv(MANIFEST, sep="\t").set_index("phenocode")
    phe_map, universe = icd10_universe()
    icd_to_phe = phe_map.assign(icd=phe_map["icd"].str.replace(".", "", regex=False)) \
                        .groupby("icd")["phecode"].apply(set).to_dict()
    phe_name = dict(zip(phe_map["phecode"], phe_map["phecode_string"]))

    rows, members, rules = [], [], {}
    for fam, names in FAMILIES.items():
        for n in names:
            rule = resolve(n, table, universe)
            rules[n] = rule
            all_icd = set().union(*rule["icd10"].values())
            phecodes = set().union(*(icd_to_phe.get(c, set()) for c in all_icd)) if all_icd else set()
            mapped = {c for c in all_icd if c in icd_to_phe}
            r = table.loc[n]
            rows.append({
                "family": fam, "endpoint": n, "longname": r["LONGNAME"],
                "cases": int(man.loc[n, "num_cases"]), "controls": int(man.loc[n, "num_controls"]),
                "n_icd10": len(all_icd), "icd10_codes": " ".join(sorted(all_icd)),
                "sources": " ".join(sorted(rule["sources"])), "mode_rule": rule["mode"],
                "atc": " ".join(sorted(rule["atc"])), "includes": " ".join(rule["includes"]),
                "case_conditions": " ; ".join(rule["conditions"]),
                "control_exclude": r["CONTROL_EXCLUDE"] if present(r["CONTROL_EXCLUDE"]) else "",
                "control_conditions": r["CONTROL_CONDITIONS"] if present(r["CONTROL_CONDITIONS"]) else "",
                "icd9_used": present(r["HD_ICD_9"]), "icd8_used": present(r["HD_ICD_8"]),
                "phecodes": " ".join(sorted(phecodes)),
                "phecode_names": " ; ".join(sorted(phe_name.get(p, p) for p in phecodes)),
                "icd10_mapped_frac": len(mapped) / len(all_icd) if all_icd else np.nan,
            })
            for src, codes in rule["icd10"].items():
                for c in codes:
                    members.append({"family": fam, "endpoint": n, "source": src, "icd10": c})
    anatomy = pd.DataFrame(rows)
    anatomy.to_csv(OUT / "endpoint_anatomy.tsv", sep="\t", index=False)
    pd.DataFrame(members).to_csv(OUT / "code_membership.tsv", sep="\t", index=False)

    pairs = []
    a = anatomy.set_index("endpoint")
    for fam, names in FAMILIES.items():
        for x, y in combinations(names, 2):
            cx = set(a.loc[x, "icd10_codes"].split()); cy = set(a.loc[y, "icd10_codes"].split())
            sx = set(a.loc[x, "sources"].split()); sy = set(a.loc[y, "sources"].split())
            px = set(a.loc[x, "phecodes"].split()); py = set(a.loc[y, "phecodes"].split())
            ctrl_same = (a.loc[x, "control_exclude"] == a.loc[y, "control_exclude"]) and \
                        (a.loc[x, "control_conditions"] == a.loc[y, "control_conditions"])
            pairs.append({
                "family": fam, "def1": x, "def2": y,
                "icd10_jaccard": jaccard(cx, cy), "source_jaccard": jaccard(sx, sy),
                "phecode_jaccard": jaccard(px, py),
                "mode_differs": a.loc[x, "mode_rule"] != a.loc[y, "mode_rule"],
                "conditions_differ": a.loc[x, "case_conditions"] != a.loc[y, "case_conditions"],
                "controls_differ": not ctrl_same,
                "log_case_ratio": abs(np.log(a.loc[x, "cases"] / a.loc[y, "cases"])),
            })
    pd.DataFrame(pairs).to_csv(OUT / "pairwise_distance.tsv", sep="\t", index=False)

    pd.set_option("display.width", 250)
    print(anatomy[["family", "endpoint", "cases", "n_icd10", "sources", "mode_rule", "atc",
                   "control_exclude", "icd10_mapped_frac", "phecode_names"]].to_string(max_colwidth=60))
    print()
    print(pd.DataFrame(pairs).round(3).to_string())


if __name__ == "__main__":
    main()
