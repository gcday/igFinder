import pandas as pd

configfile: "config.yaml"

samples = pd.read_table(config["samples"]).set_index("sample", drop=False)

rule def_all:
    input:
        config["stats_file"]
        # expand("filtered_results/mixcr/vdj_seqs/{sample}.vdj.fa", sample=samples["sample"]) #, workflow=workflow.basedir)
        
singularity: "igfinder-1.0.simg"

include: "rules/samtools_filtering.smk"
include: "rules/mixcr.smk"

def sample_to_clones():
    return {"filtered_data/mixcr/clone_summary/{0}_clone_summary.txt".format(s) : s for s in samples["sample"]}

rule gather_output:
    input:
        sample_to_clones()
        # samples = sample_to_clones()
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
    input:
        config["clones_file"]
    output:
        config["stats_file"]    
    script:
        "scripts/calc_clonality.R"

# rule report:
#     input:
#         clones=expand("results/mixcr/vdj_seqs/{sample}.vdj.fa", sample=SAMPLES)
#     output:
#         "report.html"
#     run:
#         from snakemake.utils import report
#         n_clones = 0;
#         for clones in map(open, input.clones):
#             n_clones += sum(1 for line in clones) - 1
#         report("""
#         An example variant calling workflow
#         ===================================

#         Reads were mapped to the Yeast
#         reference genome and variants were called jointly with
#         SAMtools/BCFtools.

#         This resulted in {n_clones} variants (see Table T1_).
#         """, output[0], T1=input[0])


