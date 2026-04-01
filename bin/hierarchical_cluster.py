#!/usr/bin/env python3
"""
Hierarchical clustering of pairwise structural similarity results
with automatic determination of the optimal number of clusters.

Input: TSV file from usalign_all_vs_all.bash (with header starting with #).
       Each PDB pair appears once with the best TM-score across chain combos.

Algorithm:
  1. Parse pairwise TM-scores and build a symmetric distance matrix (1 - TM)
  2. Auto mode: find the largest gap in the dendrogram merge heights
  3. Cutoff mode: cut the dendrogram at distance = 1 - cutoff
  4. Compute centroid per cluster: member with highest avg TM to cluster-mates
  5. Output in wide format (one line per cluster, centroid first)

Usage:
  # Auto-detect clusters from dendrogram gaps
  python hierarchical_cluster.py input.tsv -o clusters.dat

  # Cut at a specific TM-score threshold
  python hierarchical_cluster.py input.tsv --cutoff 0.5 -o clusters.dat

  # Force a specific number of clusters
  python hierarchical_cluster.py input.tsv --k 5 -o clusters.dat
"""

import argparse
import os
import subprocess
import sys
import time
from collections import defaultdict

import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))


def parse_args():
    p = argparse.ArgumentParser(
        description="Hierarchical clustering with automatic cluster count selection"
    )
    p.add_argument("input", help="Pairwise results TSV (from usalign_all_vs_all.bash)")
    p.add_argument(
        "--qcol", type=int, default=0,
        help="0-indexed column for query name (default: 0)"
    )
    p.add_argument(
        "--tcol", type=int, default=1,
        help="0-indexed column for target name (default: 1)"
    )
    p.add_argument(
        "--col", type=int, default=2,
        help="0-indexed column for TM-score (default: 2)"
    )
    p.add_argument(
        "--cutoff", type=float, default=None,
        help="TM-score cutoff: structures with TM >= cutoff are in the same cluster. "
             "Translates to cutting the dendrogram at height = 1 - cutoff."
    )
    p.add_argument(
        "--k", type=int, default=None,
        help="Force a specific number of clusters (skip auto-detection)"
    )
    p.add_argument(
        "--linkage", choices=["average", "complete", "single"],
        default="average",
        help="Linkage method (default: average). "
             "'ward' is not available because it requires euclidean distances."
    )
    p.add_argument(
        "-o", "--output", default=None,
        help="Output file (default: stdout)"
    )
    p.add_argument(
        "-sc", "--sort-clusters", default=None,
        help="Directory for aligned cluster PDBs"
    )
    p.add_argument(
        "--pdb-dir", default=None,
        help="PDB directory (required with --sort-clusters to resolve basenames)"
    )
    args = p.parse_args()

    if args.cutoff is not None and args.k is not None:
        p.error("--cutoff and --k are mutually exclusive")
    if args.sort_clusters and not args.pdb_dir:
        p.error("--sort-clusters requires --pdb-dir")

    return args


def load_distance_matrix(args):
    """
    Read pairwise TM-scores and build a symmetric distance matrix.
    Distance = 1 - TM_score.
    """
    col = args.col
    qcol = args.qcol
    tcol = args.tcol

    name_to_id = {}
    names = []
    scores = {}

    t0 = time.time()
    n_lines = 0

    with open(args.input, buffering=1 << 22) as fh:
        for line in fh:
            if not line or line[0] == '#':
                continue
            parts = line.rstrip('\n').split('\t')
            try:
                val = float(parts[col])
            except (IndexError, ValueError):
                continue

            n_lines += 1
            query = parts[qcol]
            target = parts[tcol]

            if query not in name_to_id:
                name_to_id[query] = len(names)
                names.append(query)
            if target not in name_to_id:
                name_to_id[target] = len(names)
                names.append(target)

            qid = name_to_id[query]
            tid = name_to_id[target]
            scores[(qid, tid)] = val
            scores[(tid, qid)] = val

    n = len(names)
    dist = np.ones((n, n), dtype=np.float64)
    np.fill_diagonal(dist, 0.0)

    for (i, j), tm in scores.items():
        dist[i, j] = 1.0 - tm

    elapsed = time.time() - t0
    print(
        f"Loaded {n_lines:,} lines in {elapsed:.1f}s, "
        f"{n:,} structures, {n_lines:,} pairs",
        file=sys.stderr,
    )
    return dist, names


def find_auto_cut(Z, n):
    """
    Find the optimal dendrogram cut by detecting the largest relative
    gap in the sorted merge heights.

    Returns the number of clusters k.
    """
    heights = Z[:, 2]

    if len(heights) < 2:
        return 2

    # Compute gaps between consecutive merge heights
    gaps = np.diff(heights)

    # Use relative gap: gap / height to handle different scales
    rel_gaps = gaps / np.maximum(heights[:-1], 1e-10)

    # Find the largest relative gap
    best_idx = np.argmax(rel_gaps)
    # Number of clusters = n - (merge_step + 1)
    k = n - best_idx - 1

    # Ensure at least 2 clusters
    k = max(k, 2)

    cut_height = heights[best_idx]
    gap_size = gaps[best_idx]
    print(
        f"Auto-cut: largest gap at merge height {cut_height:.4f} "
        f"(gap={gap_size:.4f}, relative={rel_gaps[best_idx]:.4f}), "
        f"k={k}",
        file=sys.stderr,
    )

    # Show top 5 gaps for context
    top_indices = np.argsort(rel_gaps)[::-1][:5]
    print("Top 5 dendrogram gaps:", file=sys.stderr)
    for rank, idx in enumerate(top_indices, 1):
        k_at = n - idx - 1
        marker = " <-- selected" if idx == best_idx else ""
        print(
            f"  #{rank}: height={heights[idx]:.4f}, gap={gaps[idx]:.4f}, "
            f"relative={rel_gaps[idx]:.4f}, k={k_at}{marker}",
            file=sys.stderr,
        )

    return k


def cluster_and_assign(dist, names, args):
    """Run hierarchical clustering and return clusters with centroids."""
    n = len(names)

    condensed = squareform(dist)
    Z = linkage(condensed, method=args.linkage)

    if args.k is not None:
        k = args.k
        print(f"Using user-specified k={k}", file=sys.stderr)
        labels = fcluster(Z, t=k, criterion='maxclust')
    elif args.cutoff is not None:
        cut_height = 1.0 - args.cutoff
        print(
            f"Cutting dendrogram at height={cut_height:.4f} "
            f"(TM-score cutoff={args.cutoff})",
            file=sys.stderr,
        )
        labels = fcluster(Z, t=cut_height, criterion='distance')
    else:
        k = find_auto_cut(Z, n)
        print(f"Auto-selected k={k}", file=sys.stderr)
        labels = fcluster(Z, t=k, criterion='maxclust')

    # Group members by cluster label
    clusters = defaultdict(list)
    for idx, label in enumerate(labels):
        clusters[label].append(idx)

    # Compute centroid for each cluster (member with highest avg TM to others)
    tm = 1.0 - dist
    result = []
    for label in sorted(clusters.keys()):
        members = clusters[label]
        if len(members) == 1:
            centroid = members[0]
        else:
            best_avg = -1.0
            centroid = members[0]
            for m in members:
                others = [o for o in members if o != m]
                avg_tm = np.mean([tm[m, o] for o in others])
                if avg_tm > best_avg:
                    best_avg = avg_tm
                    centroid = m
        # Put centroid first
        ordered = [centroid] + sorted(
            [m for m in members if m != centroid],
            key=lambda x: names[x]
        )
        result.append(ordered)

    # Sort clusters by size (largest first)
    result.sort(key=len, reverse=True)
    return result


def write_output(clusters, names, out):
    """Write clusters in wide format (one line per cluster, centroid first)."""
    for members in clusters:
        out.write(" ".join(names[m] for m in members) + "\n")


def main():
    args = parse_args()
    t_start = time.time()

    dist, names = load_distance_matrix(args)

    if len(names) < 2:
        print("ERROR: need at least 2 structures to cluster", file=sys.stderr)
        sys.exit(1)

    clusters = cluster_and_assign(dist, names, args)

    sizes = [len(c) for c in clusters]
    n_clusters = len(clusters)
    n_singletons = sizes.count(1)
    total = sum(sizes)
    print(
        f"Clusters: {n_clusters}, largest: {max(sizes)}, "
        f"singletons: {n_singletons}, "
        f"mean size: {total/n_clusters:.1f}, "
        f"total structures: {total}",
        file=sys.stderr,
    )

    fh = open(args.output, "w") if args.output else sys.stdout
    write_output(clusters, names, fh)
    if args.output:
        fh.close()
        print(f"Wrote {args.output}", file=sys.stderr)

    print(f"Total time: {time.time()-t_start:.1f}s", file=sys.stderr)

    if args.sort_clusters:
        if not args.output:
            print("ERROR: --sort-clusters requires --output", file=sys.stderr)
            sys.exit(1)
        align_script = os.path.join(SCRIPT_DIR, "align_clusters.bash")
        ret = subprocess.run(
            ["bash", align_script, args.output, args.sort_clusters, args.pdb_dir],
        )
        if ret.returncode != 0:
            print(
                f"ERROR: align_clusters.bash exited with code {ret.returncode}",
                file=sys.stderr,
            )
            sys.exit(ret.returncode)


if __name__ == "__main__":
    main()
