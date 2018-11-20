#!/bin/bash
#SBATCH --output=logs/master_snakemake.output
#SBATCH --error=logs/master_snakemake.output
#SBATCH --time=47:00:00

source activate snakemake

# cd /farmshare/user_data/gday/mayo/projects/igFinder
#--use-singularity --singularity-args "-B /farmshare:/farmshare" 
srun snakemake  \
--use-singularity --singularity-args "-B /farmshare:/farmshare" \
--latency-wait 100 -k  --printshellcmds -j 40 --rerun-incomplete --use-conda \
--cluster-config cluster.json  --cluster \
"sbatch --mem {cluster.mem} -c {cluster.c} --job-name {cluster.job-name} -t {cluster.time} --output .snakemake/slurm-%j.out --error .snakemake/slurm-%j.out" \
"$@"

# "sbatch  --job-name {cluster.job-name} --mem-per-cpu {cluster.mem-per-cpu} -c {cluster.c}  -t {cluster.time} "

# "sbatch  --job-name {cluster.job-name} --mem-per-cpu {cluster.mem-per-cpu} -c {cluster.c}  -t {cluster.time} "
# "slurm-%j.out"
# "sbatch  --job-name {cluster.job-name} --mem {cluster.mem} -c {cluster.c}  -t {cluster.time} --output .snakemake/slurm-%j.out --error .snakemake/slurm-%j.out"

# "sbatch  --job-name {cluster.job-name} --mem-per-cpu {cluster.mem-per-cpu} -c {cluster.c}  -t {cluster.time} --output {cluster.output} --error {cluster.error}"
