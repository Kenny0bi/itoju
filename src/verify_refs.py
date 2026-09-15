"""Verify every reference in a BibTeX file against Crossref, and every citation against the file.

For each entry with a DOI, query https://api.crossref.org/works/<doi> and compare:
  title (normalized similarity >= 0.90), first author family name, year, journal/container.
Mismatches fail loudly (nonzero exit). Entries without a DOI are listed for manual checking.
If a .tex file is given, every \\cite key must exist in the .bib and every .bib entry must be cited.

Usage: python src/verify_refs.py paper/references.bib [paper/itoju.tex]
"""
import difflib
import json
import re
import sys
import time
import unicodedata
import urllib.parse
import urllib.error
import urllib.request
from pathlib import Path

MAILTO = "itoju-reference-check"


def norm(s):
    s = unicodedata.normalize("NFKD", s or "")
    s = "".join(c for c in s if not unicodedata.combining(c))
    s = re.sub(r"\\[a-zA-Z]+|[{}\\$]", "", s)
    return re.sub(r"[^a-z0-9 ]+", " ", s.lower()).split()


def similar(a, b):
    return difflib.SequenceMatcher(None, " ".join(norm(a)), " ".join(norm(b))).ratio()


def parse_bib(text):
    entries = {}
    for m in re.finditer(r"@(\w+)\s*\{\s*([^,]+),(.*?)\n\}", text, flags=re.S):
        kind, key, body = m.group(1).lower(), m.group(2).strip(), m.group(3)
        fields = {}
        for fm in re.finditer(r"(\w+)\s*=\s*(\{(?:[^{}]|\{(?:[^{}]|\{[^{}]*\})*\})*\}|\"[^\"]*\"|\d+)", body, flags=re.S):
            v = fm.group(2).strip()
            if v[0] in "{\"":
                v = v[1:-1]
            fields[fm.group(1).lower()] = re.sub(r"\s+", " ", v).strip()
        entries[key] = {"type": kind, **fields}
    return entries


def crossref(doi):
    """Crossref record, or a DataCite record mapped to the same fields (repository DOIs such as
    Stanford Digital Repository preprints are registered with DataCite, not Crossref)."""
    url = "https://api.crossref.org/works/" + urllib.parse.quote(doi)
    req = urllib.request.Request(url, headers={"User-Agent": f"verify_refs/1.0 ({MAILTO})"})
    try:
        with urllib.request.urlopen(req, timeout=30) as r:
            return json.load(r)["message"]
    except urllib.error.HTTPError as ex:
        if ex.code != 404:
            raise
    with urllib.request.urlopen("https://api.datacite.org/dois/" + urllib.parse.quote(doi), timeout=30) as r:
        a = json.load(r)["data"]["attributes"]
    fam = lambda c: c.get("familyName") or c["name"].split(",")[0]
    return {"title": [a["titles"][0]["title"]], "author": [{"family": fam(c), "sequence": "first" if i == 0 else "additional"}
                                                             for i, c in enumerate(a["creators"])],
            "issued": {"date-parts": [[a["publicationYear"]]]}, "container-title": [a.get("publisher", "")]}


def first_family(authors_field):
    first = authors_field.split(" and ")[0].strip()
    fam = first.split(",")[0] if "," in first else first.split()[-1]
    return " ".join(norm(fam))


def check(key, e):
    doi = e.get("doi")
    if not doi:
        return "NO-DOI", "no DOI; check by hand"
    try:
        m = crossref(doi)
    except Exception as ex:  # network or 404
        return "FAIL", f"Crossref lookup failed: {ex}"
    problems = []
    title = (m.get("title") or [""])[0]
    if similar(e.get("title", ""), title) < 0.90:
        problems.append(f"title: bib '{e.get('title','')[:70]}' vs crossref '{title[:70]}'")
    year = None
    for k in ("published-print", "published-online", "issued", "created"):
        if m.get(k, {}).get("date-parts"):
            year = m[k]["date-parts"][0][0]
            break
    years = {str(m[k]["date-parts"][0][0]) for k in ("published-print", "published-online", "issued")
             if m.get(k, {}).get("date-parts") and m[k]["date-parts"][0][0]}
    if e.get("year") and e["year"] not in years:
        problems.append(f"year: bib {e['year']} vs crossref {sorted(years)}")
    authors = m.get("author") or []
    # Crossref sometimes lists a consortium (a "name" with no "family") ahead of the first person,
    # as for Grove et al. 2019. Compare against the first named person, and also accept a match with
    # the listed first author, so a consortium-first record neither hides nor fakes a mismatch.
    people = [a["family"] for a in authors if a.get("family")]
    listed_first = next((a.get("family") or a.get("name") for a in authors if a.get("sequence") == "first"),
                        (authors[0].get("family") or authors[0].get("name")) if authors else "")
    cr_first = people[0] if people else listed_first
    if e.get("author") and (cr_first or listed_first):
        bf = first_family(e["author"])
        cands = {" ".join(norm(x)) for x in (cr_first, listed_first) if x}
        if not any(bf == c or bf in c or c in bf for c in cands):
            problems.append(f"first author: bib '{bf}' vs crossref {sorted(cands)}")
    container = " ".join(m.get("container-title") or [])
    bib_container = e.get("journal") or e.get("booktitle") or ""
    if bib_container and container and similar(bib_container, container) < 0.6 \
            and not set(norm(bib_container)) <= set(norm(container)):
        problems.append(f"journal: bib '{bib_container}' vs crossref '{container}'")
    time.sleep(0.2)
    return ("FAIL", "; ".join(problems)) if problems else ("OK", f"{cr_first} {sorted(years)} {container[:40]}")


def main():
    bib = Path(sys.argv[1])
    entries = parse_bib(bib.read_text())
    failures = 0
    for key, e in entries.items():
        status, msg = check(key, e)
        failures += status == "FAIL"
        print(f"{status:6s} {key:28s} {msg}")
    if len(sys.argv) > 2:
        tex = Path(sys.argv[2]).read_text()
        cited = {k.strip() for grp in re.findall(r"\\cite[a-z]*\{([^}]*)\}", tex) for k in grp.split(",")}
        missing = cited - entries.keys()
        uncited = entries.keys() - cited
        for k in sorted(missing):
            print(f"FAIL   {k:28s} cited in tex but not in bib")
        for k in sorted(uncited):
            print(f"FAIL   {k:28s} in bib but never cited")
        failures += len(missing) + len(uncited)
    print(f"\n{len(entries)} entries, {failures} failures")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()
