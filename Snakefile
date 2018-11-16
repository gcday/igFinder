import pandas as pd

from snakemake.utils import validate, min_version
min_version("5.1.5")

configfile: "config.yaml"
validate(config, schema="schemas/config.schema.yaml")

samples = pd.read_table(config["samples"], sep = ",").set_index("sample", drop=False)

sra_table = pd.read_table(config["sra_table"], sep = "\t").set_index("Sample_Name", drop = False)
# samples = pd.read_table(config["samples"]).set_index("sample", drop=False)
# validate(samples, schema="schemas/samples.schema.yaml")

(samples_root, samples_ext) = os.path.splitext(os.path.basename(config["samples"]))
config["clones_file"] = "{}_clones.tsv".format(samples_root)
config["stats_file"] = "{}_stats.tsv".format(samples_root)

config["igblast_file"] = "{}_igblast.tsv".format(samples_root)
 
# samples_of_int = samples.query('MM_Status == ["SMM"]')
# samples_of_int = samples.query('MM_Status == ["SMM", "Untreated MM", "MGUS", "Post induction MM", "Relapsed MM", "Post ASCT"]')
sra_table = sra_table.query('Assay_Type == "RNA-Seq"')



# sample_list.remove("7216801_20110412_BM")
# samples_of_int = samples.query('')
# include: "rules/samtools_filtering.smk"
include: "rules/mixcr.smk"
include: "rules/sra.smk"

samples_of_interest = []
for sample in samples["sample"]:
    if (len(sra_table[sra_table['Sample_Name'].str.contains(str(sample))]) >= 1):
      samples_of_interest += [sample]
print(len(samples_of_interest))
samples_of_int = samples[[sample in samples_of_interest for sample in samples["sample"]]]

# samples_of_int = samples
sample_list = samples_of_int["sample"]
rule def_all:
  input:
    config["igblast_file"],
    config["clones_file"]


