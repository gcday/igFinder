#!/bin/bash
#SBATCH --output=logs/master_snakemake.output
#SBATCH --error=logs/master_snakemake.output
#SBATCH --time=2-0:0

source activate snakemake

cd /farmshare/user_data/gday/mayo/projects/igFinder

snakemake --use-singularity --singularity-args "-B /farmshare:/farmshare " --latency-wait 100 -k --printshellcmds -j 999 --rerun-incomplete --cluster-config cluster.json --cluster "sbatch --job-name {cluster.job-name} --mem-per-cpu {cluster.mem-per-cpu} -c {cluster.c}  -t {cluster.time} --output {cluster.output} --error {cluster.error}"
