#!/bin/bash
set -euo pipefail

# Align PDB structures within each cluster produced by setcover_cluster.py
#
# Input format (wide, from setcover_cluster.py --fmt wide):
#   Each line = one cluster, members space-separated
#   Each member = path/to/file.pdb:chainA:chainB
#
# Usage:
#   align_clusters.bash <clusters.tsv> <output_dir> [min_members]
#
# PDB paths in clusters.tsv are used directly (relative or absolute).
# USalign superimposes each member onto the cluster representative (first entry).

CLUSTERS="${1:?Usage: $0 <clusters.tsv> <output_dir> [min_members]}"
OUTDIR="${2:?Provide output directory}"
MIN_MEMBERS="${3:-2}"



if [ -d "$OUTDIR" ]; then
  rm -rf "$OUTDIR"
fi

mkdir -p "$OUTDIR"

# Parse entry like  path/to/file.pdb:A:B  into pdb path + chain flag
#   Sets: _pdb_path, _chain_flag, _tag
parse_entry() {
  local entry="$1"
  if [[ "$entry" == *.pdb:* ]]; then
    _pdb_path="${entry%.pdb:*}.pdb"
    local raw_chains="${entry#"$_pdb_path":}"
    # Convert colon-separated chains to comma-separated for USalign
    # Empty chain IDs (::) become _ for USalign
    local usalign_chains=""
    IFS=':' read -r -a chain_arr <<< "$raw_chains"
    for c in "${chain_arr[@]}"; do
      [[ -n "$usalign_chains" ]] && usalign_chains+=","
      if [[ -z "$c" ]]; then
        usalign_chains+="_"
      else
        usalign_chains+="$c"
      fi
    done
    _chain_flag="$usalign_chains"
    # Tag for unique output filenames: replace colons and strip trailing
    _tag="${raw_chains//[:]/}"
  else
    _pdb_path="$entry"
    _chain_flag=""
    _tag=""
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

  parse_entry "${members[0]}"
  ref_pdb="$_pdb_path"
  ref_chain_flag="$_chain_flag"
  ref_tag="$_tag"

  if [[ ! -f "$ref_pdb" ]]; then
    echo "WARNING: reference PDB not found: $ref_pdb — skipping cluster $clust_idx" >&2
    continue
  fi

  ref_base="$(basename "$ref_pdb" .pdb)"
  cp "$ref_pdb" "$cdir/${ref_base}_${ref_tag}_ref.pdb"

  # Build -chain1 flag for reference
  chain1_args=()
  [[ -n "$ref_chain_flag" ]] && chain1_args=(-chain1 "$ref_chain_flag")

  for ((i=1; i<${#members[@]}; i++)); do
    parse_entry "${members[$i]}"
    target_pdb="$_pdb_path"
    target_chain_flag="$_chain_flag"
    target_tag="$_tag"

    if [[ ! -f "$target_pdb" ]]; then
      echo "WARNING: target PDB not found: $target_pdb — skipping" >&2
      continue
    fi

    target_base="$(basename "$target_pdb" .pdb)"
    out_prefix="$cdir/${target_base}_${target_tag}_sup"

    # Build -chain2 flag for target
    chain2_args=()
    [[ -n "$target_chain_flag" ]] && chain2_args=(-chain2 "$target_chain_flag")

    USalign "$target_pdb" "$ref_pdb"  -mm 1 -ter 1 \
            "${chain1_args[@]}" "${chain2_args[@]}" \
            -o "$out_prefix" >/dev/null || \
      echo "WARNING: USalign failed for ${members[$i]}" >&2
  done

  rm -f "$cdir"/*.pml

  echo "clust$(printf '%04d' $clust_idx): ${#members[@]} members, rep=${ref_base}:${ref_tag}" >&2
done < "$CLUSTERS"

echo "Done: $clust_idx clusters written to $OUTDIR" >&2
