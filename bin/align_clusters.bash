#!/bin/bash
set -euo pipefail

# Align PDB structures within each cluster.
#
# Input format (wide):
#   Each line = one cluster, members space-separated
#   Each member = bare basename (e.g. run004_final)
#
# Usage:
#   align_clusters.bash <clusters.dat> <output_dir> <pdb_dir> [min_members]
#
# USalign superimposes each member onto the cluster representative (first entry).

CLUSTERS="${1:?Usage: $0 <clusters.dat> <output_dir> <pdb_dir> [min_members]}"
OUTDIR="${2:?Provide output directory}"
PDBDIR="${3:?Provide PDB directory}"
PDBDIR="${PDBDIR%/}"
MIN_MEMBERS="${4:-2}"

if [ -d "$OUTDIR" ]; then
  rm -rf "$OUTDIR"
fi

mkdir -p "$OUTDIR"

# Resolve bare basename to PDB file path
resolve_pdb() {
  local name="$1"
  local pdb="$PDBDIR/${name}.pdb"
  if [[ -f "$pdb" ]]; then
    echo "$pdb"
  else
    echo ""
  fi
}

clust_idx=0
while IFS= read -r line; do
  [[ -z "$line" ]] && continue

  read -r -a members <<< "$line"
  [[ ${#members[@]} -lt $MIN_MEMBERS ]] && continue

  clust_idx=$((clust_idx + 1))
  cdir="$OUTDIR/clust$(printf '%04d' $clust_idx)"
  mkdir -p "$cdir"

  ref_pdb="$(resolve_pdb "${members[0]}")"
  if [[ -z "$ref_pdb" ]]; then
    echo "WARNING: reference PDB not found: ${members[0]} — skipping cluster $clust_idx" >&2
    continue
  fi

  cp "$ref_pdb" "$cdir/${members[0]}_ref.pdb"

  for ((i=1; i<${#members[@]}; i++)); do
    target_pdb="$(resolve_pdb "${members[$i]}")"
    if [[ -z "$target_pdb" ]]; then
      echo "WARNING: target PDB not found: ${members[$i]} — skipping" >&2
      continue
    fi

    out_prefix="$cdir/${members[$i]}_sup"

    USalign "$target_pdb" "$ref_pdb" -mm 1 -ter 1 \
            -o "$out_prefix" >/dev/null || \
      echo "WARNING: USalign failed for ${members[$i]}" >&2
  done

  rm -f "$cdir"/*.pml

  echo "clust$(printf '%04d' $clust_idx): ${#members[@]} members, rep=${members[0]}" >&2
done < "$CLUSTERS"

echo "Done: $clust_idx clusters written to $OUTDIR" >&2
