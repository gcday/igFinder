#!/bin/bash
#SBATCH --output=logs/%j_master_snakemake.output
#SBATCH --error=logs/%j_master_snakemake.output
#SBATCH --time=2-0:0

source activate snakemake

cd /farmshare/user_data/gday/mayo/projects/igFinder
#--use-singularity --singularity-args "-B /farmshare:/farmshare" 
snakemake  \
--latency-wait 100 -k --printshellcmds -j 999 --rerun-incomplete --use-conda \
--cluster-config cluster.json  --cluster "sbatch  --job-name {cluster.job-name} --mem-per-cpu {cluster.mem-per-cpu} -c {cluster.c}  -t {cluster.time} --output {cluster.output} --error {cluster.error}"
