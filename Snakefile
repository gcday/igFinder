import pandas as pd
import glob
import time

from snakemake.utils import validate, min_version
min_version("5.1.5")

configfile: "config.yaml"
# validate(config, schema="schemas/config.schema.yaml")

samples = pd.read_table(config["samples"], sep = ",").set_index("sample", drop=False)
sra_table = pd.read_table(config["sra_table"], sep = "\t").set_index("Sample_Name", drop = False)
# samples = pd.read_table(config["samples"]).set_index("sample", drop=False)
# validate(samples, schema="schemas/samples.schema.yaml")

(samples_root, samples_ext) = os.path.splitext(os.path.basename(config["samples"]))
samples_root = "{}_{}".format(time.strftime("%Y%m%d"), samples_root)
config["clones_file"] = "{}_clones.tsv".format(samples_root)
config["stats_file"] = "{}_stats.tsv".format(samples_root)
config["igblast_file"] = "{}_igblast.tsv".format(samples_root)

# sample_list = [sample for sample in samples["sample"] 
											#h if (len(sra_table[sra_table['Sample_Name'].str.contains(str(sample))]) >= 1)]
# samples = samples[[sample in sample_list for sample in samples["sample"]]]
# samples_of_int = samples[[sample in samples_of_interest for sample in samples["sample"]]]
# s_to_vdjca = {s : "results/mixcr/assembled/{}.vdj.fa".format(s) for s in samples["sample"]}

# sample_list = {s for (s, path) in s_to_vdjca.items() if path in glob.glob("results/mixcr/assembled/*.vdj.fa")}
sample_list = samples["sample"][:50]
samples = samples[[sample in sample_list for sample in samples["sample"]]]

sample_list = samples["sample"]
# print(len(sample_list))



include: "rules/sra.smk"
include: "rules/mixcr.smk"
include: "rules/gather_stats.smk"
# include: "rules/samtools_filtering.smk"




 
# samples_of_int = samples.query('MM_Status == ["SMM"]')
# samples_of_int = samples.query('MM_Status == ["SMM", "Untreated MM", "MGUS", "Post induction MM", "Relapsed MM", "Post ASCT"]')
# sra_table = sra_table.query('Assay_Type == "RNA-Seq"')
# sample_list.remove("7216801_20110412_BM")
# samples_of_int = samples.query('')
# include: "rules/samtools_filtering.smk"
# s_to_vdjca = {s : "data/mixcr/aligned/{}.vdjca".format(s) for s in samples_of_int["sample"]}
# sample_list = {s for (s, path) in s_to_vdjca.items() if path in glob.glob("data/mixcr/aligned/*.vdjca")}
# s_to_clna = {s : "data/mixcr/clones/{}.clna".format(s) for s in samples_of_int["sample"]}
# sample_list = {s for (s, path) in s_to_clna.items() if path in glob.glob("data/mixcr/clones/*.clna")}
# s_to_vdjca = {s : "results/mixcr/assembled/{}.vdj.fa".format(s) for s in samples_of_int["sample"]}
# 
# sample_list = {s for (s, path) in s_to_vdjca.items() if path in glob.glob("results/mixcr/assembled/*.vdj.fa")}

# sample_list = sample_list | set(samples_of_int["sample"][:400])
# samples_of_int = samples[[sample in sample_list for sample in samples["sample"]]]

# samples_of_int = samples
# finished = "results/igblast/{sample}_igblast_output.txt"

rule def_all:
  input:
    config["igblast_file"],
    config["clones_file"]