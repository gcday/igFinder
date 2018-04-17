#!/bin/bash
#SBATCH --output=logs/master_snakemake.output
#SBATCH --error=logs/master_snakemake.output
#SBATCH --time=2-0:0

source activate snakemake

cd /farmshare/user_data/gday/mayo/projects/igFinder

snakemake --latency-wait 30 -j 999 --rerun-incomplete --cluster-config cluster.json --cluster "sbatch --job-name {cluster.job-name} --mem {cluster.mem} -n {cluster.n}  -t {cluster.time} --output {cluster.output} --error {cluster.error}"
