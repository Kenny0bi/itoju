"""Map every FinnGen endpoint definition to the OMOP standard vocabularies.

Why: FinnGen writes definitions as regular expressions over Finnish ICD-10, ICD-9 and ICD-8 codes
plus Finnish drug and reimbursement codes. To compare definitions with each other, and to rebuild
them in any OMOP CDM database, each definition is translated to standard concepts:
  - ICD-10 (WHO) codes -> "Maps to" -> SNOMED CT standard condition concepts;
  - Finnish ICD-9 codes -> ICD-9-CM (approximate: Finnish ICD-9 uses 4 digits plus an optional letter,
    e.g. 5965B; the letter is dropped and the code is dotted, 596.5) -> SNOMED;
  - ICD-8 has no OMOP vocabulary, so ICD-8 rules are counted as unmappable;
  - ATC classes (e.g. A06A, drugs for constipation) -> RxNorm ingredients via CONCEPT_ANCESTOR.

Overlap between two definitions of the same condition is then measured three ways, from coarse
to clinical: ICD-10 code Jaccard (from 02), phecodeX Jaccard (from 02), SNOMED exact-concept
Jaccard, and SNOMED hierarchy-aware Jaccard (each concept expanded to its ancestors up to
MAX_LEVELS levels, so a parent and a child concept count as related rather than disjoint).

Inputs: data/ref/omop/* (src/08a_extract_omop.sh), results/definitions/endpoint_anatomy.tsv and
pairwise_distance.tsv (src/02), the R12 endpoint file.
Outputs: results/omop/endpoint_omop.tsv, endpoint_snomed_concepts.tsv, pairwise_omop.tsv
"""
import re
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OMOP = ROOT / "data" / "ref" / "omop"
DEF = ROOT / "results" / "definitions"
ENDPOINTS = ROOT / "data" / "raw" / "endpoints" / "finngen_R12_endpoint_core_noncore_1.0.xlsx"
OUT = ROOT / "results" / "omop"
MAX_LEVELS = 2   # hierarchy-aware overlap: include ancestors within two levels

ICD9_COLS = ["HD_ICD_9", "COD_ICD_9"]
ICD8_COLS = ["HD_ICD_8", "COD_ICD_8"]


def present(v):
    return pd.notna(v) and str(v).strip() not in ("", "nan", "$!$")


def split_alternatives(pattern):
    parts, depth, cur = [], 0, ""
    for ch in pattern:
        depth += ch == "["
        depth -= ch == "]"
        if ch == "|" and depth == 0:
            parts.append(cur)
            cur = ""
        else:
            cur += ch
    parts.append(cur)
    return parts


def dot_icd10(code):
    return code if len(code) <= 3 else f"{code[:3]}.{code[3:]}"


def read(name):
    return pd.read_csv(OMOP / name, sep="\t", dtype=str, keep_default_na=False, quoting=3)


def icd9_patterns(name, table, seen=None):
    """Finnish ICD-9 patterns for an endpoint, following INCLUDE chains; plus whether ICD-8 is used."""
    seen = seen or set()
    if name in seen or name not in table.index:
        return [], False
    seen.add(name)
    r = table.loc[name]
    pats = [a.strip().lstrip("%") for c in ICD9_COLS if present(r[c]) for a in split_alternatives(str(r[c]))]
    icd8 = any(present(r[c]) for c in ICD8_COLS)
    if present(r["INCLUDE"]):
        for child in str(r["INCLUDE"]).split("|"):
            p, e8 = icd9_patterns(child.strip(), table, seen)
            pats += p
            icd8 |= e8
    return pats, icd8


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    src = read("source_concepts.tsv")
    maps = read("maps_to.tsv")
    tgt = read("target_concepts.tsv").set_index("concept_id")
    anc = read("target_ancestors.tsv")
    atc_desc = read("atc_descendants.tsv")
    related = read("related_concepts.tsv").set_index("concept_id")
    versions = read("vocabulary_versions.tsv") if (OMOP / "vocabulary_versions.tsv").exists() else None

    maps_to = maps.groupby("concept_id_1")["concept_id_2"].apply(set).to_dict()
    icd10 = src[src.vocabulary_id == "ICD10"].set_index("concept_code")["concept_id"].to_dict()
    icd9 = src[src.vocabulary_id == "ICD9CM"]
    icd9_undotted = dict(zip(icd9.concept_code.str.replace(".", "", regex=False), icd9.concept_id))
    atc = src[src.vocabulary_id == "ATC"].set_index("concept_code")["concept_id"].to_dict()

    anc = anc[anc.min_levels_of_separation.astype(int) <= MAX_LEVELS]
    closure = anc.groupby("descendant_concept_id")["ancestor_concept_id"].apply(set).to_dict()
    snomed_std = set(tgt.index[(tgt.vocabulary_id == "SNOMED") & (tgt.standard_concept == "S")])

    anatomy = pd.read_csv(DEF / "endpoint_anatomy.tsv", sep="\t")
    table = pd.read_excel(ENDPOINTS).set_index("NAME")

    rows, concept_rows, sets = [], [], {}
    for _, a in anatomy.iterrows():
        e = a.endpoint
        codes = str(a.icd10_codes).split() if pd.notna(a.icd10_codes) else []
        cids = [icd10.get(dot_icd10(c)) for c in codes]
        found = [c for c in cids if c]
        snomed = set().union(*(maps_to.get(c, set()) for c in found)) & snomed_std if found else set()

        pats, uses_icd8 = icd9_patterns(e, table)
        icd9_hits = set()
        for p in pats:
            rx = re.compile("^" + re.sub(r"[A-Z]$", "", p.replace(".", "")))
            icd9_hits |= {cid for code, cid in icd9_undotted.items() if rx.match(code)}
        snomed9 = set().union(*(maps_to.get(c, set()) for c in icd9_hits)) & snomed_std if icd9_hits else set()

        atc_codes = str(a.atc).split() if pd.notna(a.atc) and str(a.atc).strip() else []
        ingredients = set()
        for code in atc_codes:
            cid = atc.get(code)
            if cid:
                desc = set(atc_desc.loc[atc_desc.ancestor_concept_id == cid, "descendant_concept_id"])
                ingredients |= {d for d in desc if d in related.index
                                and related.loc[d, "vocabulary_id"] == "RxNorm"
                                and related.loc[d, "concept_class_id"] == "Ingredient"}

        all_snomed = snomed | snomed9
        sets[e] = all_snomed
        rows.append({
            "family": a.family, "endpoint": e, "cases": a.cases,
            "icd10_codes": len(codes), "icd10_in_omop": len(found),
            "icd10_mapped_to_snomed": sum(1 for c in found if maps_to.get(c, set()) & snomed_std),
            "snomed_from_icd10": len(snomed),
            "icd9_patterns": len(pats), "icd9cm_concepts_matched": len(icd9_hits), "snomed_from_icd9": len(snomed9),
            "snomed_total": len(all_snomed), "uses_icd8_unmappable": uses_icd8,
            "atc_classes": " ".join(atc_codes), "rxnorm_ingredients": len(ingredients),
        })
        for c in sorted(all_snomed):
            concept_rows.append({"endpoint": e, "snomed_concept_id": c, "from_icd10": c in snomed,
                                 "domain_id": tgt.loc[c, "domain_id"]})

    omop = pd.DataFrame(rows)
    omop.to_csv(OUT / "endpoint_omop.tsv", sep="\t", index=False)
    pd.DataFrame(concept_rows).to_csv(OUT / "endpoint_snomed_concepts.tsv", sep="\t", index=False)

    dist = pd.read_csv(DEF / "pairwise_distance.tsv", sep="\t")
    out = []
    for _, d in dist.iterrows():
        x, y = sets.get(d.def1, set()), sets.get(d.def2, set())
        cx = set().union(x, *(closure.get(c, set()) for c in x))
        cy = set().union(y, *(closure.get(c, set()) for c in y))
        jac = lambda p, q: len(p & q) / len(p | q) if (p | q) else np.nan
        out.append({**d.to_dict(), "snomed_jaccard": jac(x, y), "snomed_hier_jaccard": jac(cx, cy),
                    "snomed_n1": len(x), "snomed_n2": len(y)})
    pw = pd.DataFrame(out)
    pw.to_csv(OUT / "pairwise_omop.tsv", sep="\t", index=False)

    pd.set_option("display.width", 250)
    if versions is not None:
        print(versions.to_string(index=False))
    print(omop.to_string(index=False))
    print(pw[["family", "def1", "def2", "icd10_jaccard", "phecode_jaccard", "snomed_jaccard",
              "snomed_hier_jaccard"]].round(3).to_string(index=False))


if __name__ == "__main__":
    main()
