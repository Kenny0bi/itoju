#!/bin/bash
# Pull only the slices of the OMOP vocabularies (Athena download) that itoju needs.
#
# The full download is about 4.8 GB and lives outside the repo. The slices go to data/ref/omop/,
# which is git-excluded: SNOMED CT content is licensed and must not be redistributed, so only
# concept IDs and counts derived from it are published.
#
# Usage: bash src/08a_extract_omop.sh /path/to/vocabulary_download_folder
set -euo pipefail
ATHENA="$1"
cd "$(dirname "$0")/.."
OUT=data/ref/omop
mkdir -p "$OUT"

# 1. Source concepts: WHO ICD-10, ICD-9-CM, ATC.
awk -F'\t' 'NR==1 || $4=="ICD10" || $4=="ICD9CM" || $4=="ATC"' "$ATHENA/CONCEPT.csv" > "$OUT/source_concepts.tsv"

# 2. "Maps to" relationships from those source concepts to standard concepts (valid rows only).
awk -F'\t' 'NR==FNR{if(FNR>1) s[$1]=1; next} FNR==1 || (($1 in s) && $3=="Maps to" && $6=="")' \
    "$OUT/source_concepts.tsv" "$ATHENA/CONCEPT_RELATIONSHIP.csv" > "$OUT/maps_to.tsv"

# 3. The standard concepts those map to.
awk -F'\t' 'NR==FNR{if(FNR>1) t[$2]=1; next} FNR==1 || ($1 in t)' \
    "$OUT/maps_to.tsv" "$ATHENA/CONCEPT.csv" > "$OUT/target_concepts.tsv"

# 4. Ancestors of every mapped standard concept (for hierarchy-aware overlap), and the drug
#    descendants of the ATC classes (to count the RxNorm ingredients a drug rule captures).
awk -F'\t' 'NR==FNR{if(FNR>1) t[$2]=1; next} FNR==1 || ($2 in t)' \
    "$OUT/maps_to.tsv" "$ATHENA/CONCEPT_ANCESTOR.csv" > "$OUT/target_ancestors.tsv"
awk -F'\t' 'NR==FNR{if(FNR>1 && $4=="ATC") a[$1]=1; next} FNR==1 || ($1 in a)' \
    "$OUT/source_concepts.tsv" "$ATHENA/CONCEPT_ANCESTOR.csv" > "$OUT/atc_descendants.tsv"

# 5. Names and classes for every concept referenced by the ancestor and descendant slices.
awk -F'\t' 'FILENAME==ARGV[1]{if(FNR>1) k[$1]=1; next} FILENAME==ARGV[2]{if(FNR>1) k[$2]=1; next} FNR==1 || ($1 in k)' \
    "$OUT/target_ancestors.tsv" "$OUT/atc_descendants.tsv" "$ATHENA/CONCEPT.csv" > "$OUT/related_concepts.tsv"

grep -E "^(SNOMED|ICD10|ICD9CM|ATC|RxNorm)	" "$ATHENA/VOCABULARY.csv" | cut -f1,4 > "$OUT/vocabulary_versions.tsv"
wc -l "$OUT"/*.tsv
touch "$OUT/extract.complete"
