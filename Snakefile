import pandas as pd
# SHORT_SAMPLES = ["s_27779"]

# SAMPLES = ["s_27779", "s_27780", "s_27817", "s_27818", "s_28016", "s_28017", "s_28019", "s_28240", "s_28242", "s_28243", "s_28246",
#            "s_28248", "s_28249", "s_28251", "s_28253", "s_28255", "s_28257", "s_28259", "s_28261", "s_28262", "s_28264", "s_28459", 
#            "s_28460", "s_28461", "s_28462", "s_28465", "s_28466", "s_28471", "s_28472", "s_28473", "s_28803", "s_28804", "s_28805", 
#            "s_28806", "s_28807", "s_28808", "s_28809", "s_28810", "s_28811", "s_28812", "s_28813", "s_28814", "s_28815", "s_28816", 
#            "s_28817", "s_28818", "s_28819", "s_28827", "s_28828", "s_28828r", "s_28829r", "s_28830", "s_29019", "s_29020", "s_29021", 
#            "s_29022", "s_29023", "s_29024", "s_29025", "s_29026", "s_29027", "s_29028", "s_29029", "s_29030", "s_29031", "s_29032", 
#            "s_29033", "s_29034", "s_29035", "s_29036", "s_29037", "s_29038", "s_29039", "s_29040", "s_29041", "s_29042", "s_29043", 
#            "s_29044", "s_29045", "s_29046", "s_29047", "s_29048", "s_29049", "s_29050", "s_29106", "s_29107", "s_29108", "s_29109", 
#            "s_29110", "s_29111", "s_29112", "s_29113", "s_29114", "s_29115", "s_29165", "s_29166", "s_29167", "s_29168", "s_29169", 
#            "s_29170", "s_29171", "s_29172", "s_29173", "s_29174", "s_29175", "s_29176", "s_29177", "s_29178", "s_29179", "s_29180", 
#            "s_29181", "s_29184", "s_29187", "s_29188", "s_29189", "s_29190", "s_29191", "s_29192", "s_29194", "s_29195", "s_29196", 
#            "s_29198", "s_29199", "s_29200", "s_29201", "s_29202", "s_29203", "s_29204", "s_29208", "s_29210", "s_29212", "s_29215", 
#            "s_29216", "s_29217", "s_29218", "s_29219", "s_29227", "s_29228", "s_GM15510-1", "s_GM15510-2", "s_GM15510-3", "s_GM15510-4", 
#            "s_GM19240-1", "s_GM19240-2", "s_GM19240-3", "s_MC1286PE1", "s_MC1286PE2", "s_MC1286PE3", "s_MC1286PE5", "s_MC1286PE7"]
# WORKFLOW = workflow.basedir

configfile: "config.yaml"

samples = pd.read_table(config["samples"]).set_index("sample", drop=False)
print(samples["sample"])

rule def_all:
    input:
        expand("filtered_results/mixcr/vdj_seqs/{sample}.vdj.fa", sample=samples["sample"]) #, workflow=workflow.basedir)


singularity: "docker://gcday/igfinder-1.0"

# singularity: "docker://gcday/igfinder:1.0"
# rule samtools_namesort_fastq:
#     input:
#         "data/bam_files/{sample}.bam"
#     output:
#         one=temp("data/fastq/{sample}_1.fastq"),
#         two=temp("data/fastq/{sample}_2.fastq"),
#         singleton=temp("data/fastq/{sample}_singleton.fastq")
#     threads: 16
#     shell:
#         "samtools sort -n -m 2000M -@ {threads} {input} " 
#         "| samtools fastq -N -@ {threads} "
#         "-1 {output.one} -2 {output.two} -0 {output.singleton} -"
        # ig_bed="Ig_exons.sorted.extended300.merged.bed"
        # bam="data/bam_files/{sample}.bam",

def get_bam(wildcards):
    return samples.loc[(wildcards.sample), ["bam"]].dropna()

rule samtools_filter:
    input:
        bam=get_bam,
        ig_bed=config["ig_bed"]
    output:
        merged=temp("data/filtered_bam/{sample}.bam"),
        temp1=temp("data/filtered_bam/temp_1{sample}.bam"),
        temp2=temp("data/filtered_bam/temp_2{sample}.bam")
    threads: 16
    shell:
        "samtools view -uh -@ {threads} -f 13 {input.bam} | samtools sort -@ {threads} -o {output.temp1} &&"
        "samtools view -uh -@ {threads} -f 1 -M -L {input.ig_bed} {input.bam} | samtools sort -@ {threads} -o {output.temp2} &&"
        "samtools merge {output.merged} {output.temp1} {output.temp2} "

rule samtools_sort:
    input:
        rules.samtools_filter.output.merged
    output:
        temp("data/filtered_namesorted_reads/{sample}.bam")
    threads: 16
    shell:
        "samtools sort -t /tmp -n -m 2000M -@ {threads} -o {output} {input}"

rule samtools_fastq:
    input:
        rules.samtools_sort.output
    output:
        one=temp("data/fastq/{sample}_1.fastq"),
        two=temp("data/fastq/{sample}_2.fastq"),
        singleton=temp("data/fastq/{sample}_singleton.fastq")
    threads: 16
    shell:
        "samtools fastq -N -@ {threads} "
        "-1 {output.one} -2 {output.two} -s {output.singleton} {input}"


rule mixcr_align:
    input:
        rules.samtools_fastq.output.one,
        rules.samtools_fastq.output.two
    output:
        "filtered_data/mixcr/aligned/{sample}.vdjca"
    threads: 16
    shell: 
        "mixcr align -t {threads} -g -a -f -p rna-seq -s hsa "
        "-OvParameters.geneFeatureToAlign=VGeneWithP "
        "-OallowPartialAlignments=true "
        "{input} {output}"

rule mixcr_assemble:
    input:
        rules.mixcr_align.output
    output:
        index="data/mixcr/index/{sample}.index",
        clones="data/mixcr/clones/{sample}.clns"
    shell: 
        "mixcr assemble -f -i {output.index} {input} {output.clones}"

rule mixcr_export:
    input:
        rules.mixcr_assemble.output.clones
    output:
        "filtered_data/mixcr/clones/{sample}.txt"
    shell:
        "mixcr exportClones -cloneId -count -fraction -vGene -dGene -jGene "
        "-aaFeature CDR3 -vBestIdentityPercent -vIdentityPercents "
        "-jBestIdentityPercent -jIdentityPercents -nFeature CDR3 "
        "-avrgFeatureQuality CDR3 -minFeatureQuality CDR3 "
        "{input} {output}"

rule mixcr_export_sig_clones:
    input:
        rules.mixcr_export.output,
        rules.mixcr_align.output,
        rules.mixcr_assemble.output.index
    output:
        "ftmp/mixcr/{sample}",
        "filtered_results/mixcr/vdj_seqs/{sample}.vdj.fa"
    threads: 16
    shell:
        "scripts/assemble_clones.sh {input} {output} {threads} {wildcards.sample}"

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


