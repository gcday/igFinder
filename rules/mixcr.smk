def get_srr(wildcards):
    return list(sra_table[sra_table['Sample_Name'].str.contains(str(wildcards.sample))]['Run'])[0]
rule download:
    params:
        srr=get_srr,
        sra_dir = config["sra_dir"],
        fq1_path = os.path.abspath("data/temp/fastq/{sample}_1.fastq"),
        fq2_path = os.path.abspath("data/temp/fastq/{sample}_2.fastq"),
        log_path = os.path.abspath("logs/download/{sample}.log")
    output:
        fastq1=temp("data/temp/fastq/{sample}_1.fastq"),
        fastq2=temp("data/temp/fastq/{sample}_2.fastq"),
        # bam_file="data/bam/{sample}.bam"
    log:
        "logs/download/{sample}.log"
    group: "igFinder"
    threads: 16
    shell:
        "cd {params.sra_dir}; fasterq-dump -p {params.srr} 1>{params.log_path} 2>&1 && "
        "mv {params.srr}_1.fastq {params.fq1_path} && mv {params.srr}_2.fastq {params.fq2_path} "

rule mixcr_align:
  input:
    rules.download.output.fastq1,
    rules.download.output.fastq2
  #conda: "../envs/igFinder.yaml"
  log:
    "logs/mixcr_align/{sample}.log"
  output:
    "data/mixcr/aligned/{sample}.vdjca"
  threads: 16
  group: "igFinder"
  shell: 
    "mixcr align -t {threads} -p kAligner2 -f -s hsa "
    "-OallowPartialAlignments=true -OsaveOriginalReads=true "
    "{input} {output} 1>{log} 2>&1"
    # "mixcr align -t {threads} -g -a -f -p rna-seq -s hsa "


rule mixcr_assemble:
  input:
    rules.mixcr_align.output
  #conda: "../envs/igFinder.yaml"
  output:
    clones="data/mixcr/clones/{sample}.clna"
  log:
    "logs/mixcr_assemble/{sample}.log"
  group: "igFinder"
  threads: 8
  shell: 
    "mixcr assemble -a -f -t {threads} {input} {output.clones} 1>{log} 2>&1"

rule mixcr_export:
  input:
    rules.mixcr_assemble.output.clones
  #conda: "../envs/igFinder.yaml"
  output:
    short="data/mixcr/clone_summary/{sample}_clone_summary.txt",
    func_clones="data/mixcr/func_clones/{sample}_func_clones.txt",
    full="data/mixcr/full_clones/{sample}_full_clones.txt"
  params:
    short = "-cloneId -count -fraction -vGene -dGene -jGene -aaFeature CDR3 -vBestIdentityPercent -vIdentityPercents -jBestIdentityPercent -jIdentityPercents -nFeature CDR3 -avrgFeatureQuality CDR3 -minFeatureQuality CDR3"
  log:
    "logs/mixcr_export/{sample}.log"
  group: "igFinder"
  shell:
    "mixcr exportClones -t -o {params.short} "
    "{input} {output.func_clones}; "
    "mixcr exportClones {params.short} "
    "{input} {output.short}; "
    "mixcr exportClones --preset full {input} {output.full}"

rule new_mixcr_export_sig:
  input:
    clones=rules.mixcr_export.output.func_clones,
    clna=rules.mixcr_assemble.output
  output:
    tempdir = directory("tmp/mixcr/{sample}"),
    contigs = "results/mixcr/top_func_seq/{sample}.vdj.fa"
  log:
    "logs/mixcr_export_sig_clones/{sample}.log"
  group: "igFinder"
  threads: 1
  script:
    "../scripts/assemble_clones.py"

rule igblast:
  input:
    contigs="results/mixcr/top_func_seq/{sample}.vdj.fa",
    db_V = config["igblast"]["germline_db_V"],
    db_J = config["igblast"]["germline_db_J"],
    db_D = config["igblast"]["germline_db_D"],
    aux_file = config["igblast"]["auxiliary_file"]
  output:
    "results/igblast/{sample}_igblast_output.txt"
  group: "igFinder"
  shell:
    "if [ -s {input.contigs} ]\n"
    "then\n"
    "    igblastn -germline_db_V {input.db_V} -germline_db_J {input.db_J} "
    "-germline_db_D {input.db_D} -query {input.contigs} -outfmt 7 "
    "-auxiliary_data {input.aux_file} > {output} \n"
    "else\n"
    "    touch {output}; exit\n"
    "fi" 