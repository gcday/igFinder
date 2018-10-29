import pandas as pd

from snakemake.utils import validate, min_version
min_version("5.1.5")

configfile: "config.yaml"
validate(config, schema="schemas/config.schema.yaml")
samples = pd.read_table(config["samples"]).set_index("sample", drop=False)
validate(samples, schema="schemas/samples.schema.yaml")

(samples_root, samples_ext) = os.path.splitext(os.path.basename(config["samples"]))
config["clones_file"] = "{}_clones.tsv".format(samples_root)
config["stats_file"] = "{}_stats.tsv".format(samples_root)

tumor_list = samples.query('MM_Status == ["SMM"]')["sample"]


include: "rules/samtools_filtering.smk"
include: "rules/mixcr.smk"

rule def_all:
  input:
      expand("results/igblast/{sample}_igblast_output.txt", sample=tumor_list)

    # expand(["data/fastq/{sample}_1.fastq", "data/fastq/{sample}_2.fastq"], sample=samples["sample"][1])

    # expand("results/igblast/{sample}_igblast_output.txt", sample=samples["sample"][1])
    # config["stats_file"],
    # expand("results/mixcr/vdj_seqs/{sample}.vdj.fa", sample=samples["sample"]) #, workflow=workflow.basedir)

    # "read_counts.tsv",


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
      summary.write("sample\tread_count\n")
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
    config["clones_file"],
    "read_counts.tsv"
  output:
    config["stats_file"]    
  script:
    "scripts/calc_clonality.R"
