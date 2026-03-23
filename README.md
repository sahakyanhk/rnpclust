# rnpclust

Wrapper for strictural clustering of ribonucleoproteins using USalign.

## Prerequisites

- [USalign](https://github.com/pylelab/USalign)
- [GNU parallel](https://www.gnu.org/software/parallel/)
- Python 3.6+

## Installation

```bash
git clone https://github.com/sahakyanhku/rnpclust
cd rnpclust
bash install.sh
```

This checks prerequisites and prints a `PATH` export line for your shell profile.

## Usage

```bash
# Full pipeline: all-vs-all alignment → clustering → structural superposition
uniclust <pdb_dir> [output_dir] [cutoff]

# Example
uniclust my_structures/ results/ 0.5
```

### Output

```
results/
  all_vs_all.tsv           # pairwise USalign scores
  clusters.dat             # cluster assignments (wide format)
  aligned_clusters/        # superimposed PDBs per cluster
    clust0001/
    clust0002/
    ...
```

## Step-by-step

Each step can be run independently:

```bash
# All-vs-all structural alignment
bin/usalign_all_vs_all.bash <pdb_dir> [output.tsv]

# Clustering only (works with any pairwise score TSV)
python bin/setcover_cluster.py input.tsv --col 2 --cutoff 0.5 -o clusters.dat

# Align within clusters
bin/align_clusters.bash clusters.dat output_dir/
```

