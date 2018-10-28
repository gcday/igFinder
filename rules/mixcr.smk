rule mixcr_align:
  input:
    rules.samtools_fastq.output.one,
    rules.samtools_fastq.output.two
  conda: "../envs/igFinder.yaml"
  log:
    "logs/mixcr_align/{sample}.log"
  output:
    "data/mixcr/aligned/{sample}.vdjca"
  threads: 16
  #group: "igFinder"
  shell: 
    "mixcr align -t {threads} -g -a -f -p rna-seq -s hsa "
    "-OallowPartialAlignments=true "
    "-OvParameters.geneFeatureToAlign=VGeneWithP "
    "{input} {output}"

rule mixcr_assemble:
  input:
    rules.mixcr_align.output
  conda: "../envs/igFinder.yaml"
  output:
    index="data/mixcr/index/{sample}.index",
    clones="data/mixcr/clones/{sample}.clns"
  log:
    "logs/mixcr_assemble/{sample}.log"
    #group: "igFinder"
  threads: 4
  shell: 
    "mixcr assemble -f -t {threads} -i {output.index} {input} {output.clones}"

rule mixcr_export:
  input:
    rules.mixcr_assemble.output.clones
  conda: "../envs/igFinder.yaml"
  output:
    short="data/mixcr/clone_summary/{sample}_clone_summary.txt",
    full="data/mixcr/full_clones/{sample}_full_clones.txt"
  log:
    "logs/mixcr_export/{sample}.log"
    #group: "igFinder"
  shell:
    "mixcr exportClones -cloneId -count -fraction -vGene -dGene -jGene "
    "-aaFeature CDR3 -vBestIdentityPercent -vIdentityPercents "
    "-jBestIdentityPercent -jIdentityPercents -nFeature CDR3 "
    "-avrgFeatureQuality CDR3 -minFeatureQuality CDR3 "
    "{input} {output.short}; "
    "mixcr exportClones --preset full {input} {output.full}"

rule mixcr_export_sig_clones:
  input:
    rules.mixcr_export.output.short,
    rules.mixcr_align.output,
    rules.mixcr_assemble.output.index
  conda: "../envs/igFinder.yaml"
  log:
    "logs/mixcr_export_sig_clones/{sample}.log"
  output:
    directory("tmp/mixcr/{sample}"),
    "results/mixcr/vdj_seqs/{sample}.vdj.fa"
  threads: 16
  #group: "igFinder"
  shell:
    "scripts/assemble_clones.sh {input} {output} {threads} "
    "{wildcards.sample}"

rule igblast:
  input:
    "results/mixcr/vdj_seqs/{sample}.vdj.fa"

  ""

# rule mixcr_join_sig_clones:
#   input:
#     rules.mixcr_export.output.short,
#     rules.mixcr_export_sig_clones.output,


# min_reads_for_assembly = 3
# min_fraction_for_assembly = 0.05
# def get_sig_clones(wildcards, input):
#   clones = pd.read_table(input["clones"]).set_index("cloneID", drop=False)
#   sig_clones = clones.loc[int(clones["cloneCount"]) >= min_reads_for_assembly |
#                     int(clones["cloneFraction"]) >= min_fraction_for_assembly]
#   return " ".join(sig_clones["cloneID"])

# rule export_reads_sig_clones:
#   input:
#     clones=rules.mixcr_export.output.short,
#     vdjca=rules.mixcr_align.output,
#     index=rules.mixcr_assemble.output.index
#   params:
#     export_name="data/mixcr_reads/{sample}/reads.fastq",
#     clone_list=get_sig_clones(wildcards, input)
#   output:
#     r1=dynamic("data/mixcr_reads/{sample}/reads_cln{cloneID}_R1.fastq"),
#     r2=dynamic("data/mixcr_reads/{sample}/reads_cln{cloneID}_R2.fastq")
#   log:
#     "logs/mixcr_export_reads/{sample}.log"
#   run:
#     clones = pd.read_table(input["clones"]).set_index("cloneID", drop=False)
#     sig_clones = clones.loc[clones["cloneCount"] >= min_reads_for_assembly |
#                     clones["cloneFraction"] >= min_fraction_for_assembly]
#     sig_list = " ".join(sig_clones["cloneID"])
#     shell:
#     "mixcr exportReadsForClones {input.index} {input.vdjca} "
#       " {params.clone_list} {params.export_name} "

#     # shell(
# rule assemble_sig_clone:
#   shadow: "shallow"
#   input:
#     r1=dynamic("data/mixcr_reads/{sample}/reads_cln{cloneID}_R1.fastq"),
#     r2=dynamic("data/mixcr_reads/{sample}/reads_cln{cloneID}_R2.fastq")
#   # output:
#   params:
#     vOpt_dir="data/vOpt/{sample}/cln_{cloneID}_vOptOut"
#   output:
#     contigs=touch("data/vOpt/{sample}/cln_{cloneID}_vOptOut/contigs.fa")
#     # vOpt_d
#   threads: 16
#   shell:
#     "VelvetOptimiser.pl -v -s 51 -e 151 -x 4 -t {threads} -c n50*tbp "
#     "-k n50*tbp -a -d {params.vOpt_dir} "
#     "-f \" -fastq -shortPaired -separate {input.r1} {input.r2} \" "

# # rule gather_contigs:
#   input:
#     dynamic("data/vOpt/{sample}/cln_{cloneID}_vOptOut/contigs.fa")
#   output:
#     "results/contigs/{sample}/cln_{cloneID}.fastq"
#   run:
    
# rule rename_contigs:
#   input:
#     dynamic("data/vOpt/{sample}/cln_{cloneID}_vOptOut/contigs.fa")
#   output:
#     "data/vOpt/{sample}/cln_{cloneID}_vOptOut/contigs.fa"
# rule mixcr_export_sig_clones_python:
#   input:
#     clones=rules.mixcr_export.output.short,
#     vdjca=rules.mixcr_align.output,
#     index=rules.mixcr_assemble.output.index
#   conda: "../envs/igFinder.yaml"
#   log:
#     "logs/mixcr_export_sig_clones/{sample}.log"
#   output:
#     "tmp/mixcr/{sample}",
#     "data/mixcr/reads/{sample}.fa"
#     "results/mixcr/vdj_seqs/{sample}.vdj.fa"
#   threads: 16
#     #group: "igFinder"
#   shell:
#   run:
    
#     for clone in sig_clones:
#       clone_id = int(clone["cloneID"])
#       shell("mixcr exportReadsForClones {input.index} {input.vdjca} {cloneID} ${OUT}/reads.fastq")
#         "scripts/assemble_clones.sh {input} {output} {threads} "
#         "{wildcards.sample}"