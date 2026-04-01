#!/bin/bash
set -euo pipefail

DIR="${1:?Usage: $0 <pdb_dir> [output_file]}"
DIR="${DIR%/}"
SUFFIX=".pdb"

TMPLIST=$(mktemp /tmp/usalign_list.XXXXXX)
trap 'rm -f "$TMPLIST"' EXIT

# List basenames without suffix
for f in "$DIR"/*"$SUFFIX"; do
  basename "$f" "$SUFFIX"
done > "$TMPLIST"

N=$(wc -l < "$TMPLIST")
echo "Found $N structures → $(( N*(N-1)/2 )) pairs" >&2

header="#PDBchain1\tPDBchain2\tTM\tTM1\tTM2\tRMSD\tID1\tID2\tIDali\tL1\tL2\tLali"

run_alignments() {
  local total=$(( N*(N-1)/2 ))
  local counter=0
  awk '{a[NR]=$1} END{for(i=1;i<NR;i++) for(j=i+1;j<=NR;j++) print a[i],a[j]}' \
      "$TMPLIST" \
    | parallel --will-cite --colsep ' ' -j "$(nproc)" \
      "USalign '${DIR}'/{1}${SUFFIX} '${DIR}'/{2}${SUFFIX} -mm 1 -ter 1 -outfmt 2 2>/dev/null | grep -v '^#'" \
    | while IFS= read -r line; do
        counter=$((counter + 1))
        printf '\r%d / %d completed' "$counter" "$total" >&2
        printf '%s\n' "$line"
      done
  printf '\n' >&2
}

# Add TM=max(TM1,TM2) column after PDBchain2
add_tm_col() {
  echo -e "$header"
  awk -F'\t' 'BEGIN{OFS="\t"} {tm=($3>$4)?$3:$4; print $1,$2,tm,$3,$4,$5,$6,$7,$8,$9,$10,$11}'
}

if [[ -n "${2:-}" ]]; then
  run_alignments | add_tm_col > "$2"
  echo "Done → $2" >&2
else
  run_alignments | add_tm_col
fi



