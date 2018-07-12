import pandas as pd

from snakemake.utils import validate, min_version
min_version("5.1.5")

configfile: "config.yaml"
validate(config, schema="schemas/config.schema.yaml")
samples = pd.read_table(config["samples"]).set_index("sample", drop=False)
validate(samples, schema="schemas/samples.schema.yaml")

rule def_all:
  input:
    config["stats_file"],
    "read_counts.tsv"
    #expand("results/mixcr/vdj_seqs/{sample}.vdj.fa", sample=samples["sample"]) #, workflow=workflow.basedir)

include: "rules/samtools_filtering.smk"
include: "rules/mixcr.smk"

def sample_to_clones():
  return {"data/mixcr/clone_summary/{0}_clone_summary.txt".format(s) : s for s in samples["sample"]}
def sample_to_read_counts():
  return {"data/read_counts/{0}.txt".format(s) : s for s in samples["sample"]}


rule gather_read_counts:
  input:
    sample_to_read_counts()
  output:
    "read_counts.tsv"
  run:
    files = sample_to_read_counts()
    with open(output[0], "w+") as summary:
      summary.write("sample\tread_count")
      for sample in files:
        with open(sample) as file:
          count = file.readline()
          summary.write(files[sample] + "\t" + count)
          # if len(clones) < 2:
            # summary.write(files[sample] + "\n")
            # continue
          # for line in clones[1:]:
            # summary.write(files[sample] + "\t" + line)

rule gather_output:
  input:
    sample_to_clones()
  output:
    config["clones_file"]
  run:
    files = sample_to_clones()
    with open(output[0], "w+") as summary:
      first = True
      for sample in files:
        with open(sample) as file:
          clones = file.readlines()
          if (first):
            summary.write("sample\t" + clones[0])
            first = False
          if len(clones) < 2:
            summary.write(files[sample] + "\n")
            continue
          for line in clones[1:]:
            summary.write(files[sample] + "\t" + line)
            
rule calc_clonality:
  conda: "envs/R.yaml"
  input:
    config["clones_file"]
  output:
    config["stats_file"]    
  script:
    "scripts/calc_clonality.R"
