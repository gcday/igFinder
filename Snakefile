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

config["igblast_file"] = "{}_igblast.tsv".format(samples_root)
 
# samples_of_int = samples.query('MM_Status == ["SMM"]')
# samples_of_int = samples.query('MM_Status == ["SMM", "Untreated MM", "MGUS", "Post induction MM", "Relapsed MM", "Post ASCT"]')


samples_of_int = samples
sample_list = samples_of_int["sample"]
# sample_list.remove("7216801_20110412_BM")
# samples_of_int = samples.query('')
include: "rules/samtools_filtering.smk"
include: "rules/mixcr.smk"

rule def_all:
  input:
    config["igblast_file"],
    config["clones_file"]


def sample_to_clones():
  return {"data/mixcr/clone_summary/{0}_clone_summary.txt".format(s) : s for s in sample_list}
def sample_to_read_counts():
  return {"data/read_counts/{0}.txt".format(s) : s for s in sample_list}

def sample_to_igblast():
  return {"results/igblast/{0}_igblast_output.txt".format(s) : s for s in sample_list}

def sample_to_igblast():
  return {s :"results/igblast/{0}_igblast_output.txt".format(s) for s in sample_list}


rule gather_igblast:
  input:
    expand("results/igblast/{sample}_igblast_output.txt", sample = sample_list)
  output:
    config["igblast_file"]
  run:
    files = sample_to_igblast()
    igblast_samples = samples_of_int.copy()
    max_V_ident = []
    for sample in samples_of_int["sample"]:
      with open(files[sample]) as igblast:
        max_ident = ""
        for line in igblast.readlines():
          if "Total\tN/A" in line:
            split = line.strip().split()
            curr_ident = float(split[-1])
            max_ident = str(curr_ident)
            if (curr_ident < 99.0):
              break
        max_V_ident += [max_ident]
    igblast_samples["max_V_ident"] = max_V_ident
    igblast_samples = igblast_samples.drop(["sample", "mutect2.vcf"], axis = 1)
    igblast_samples.to_csv(output[0], sep = "\t")

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
