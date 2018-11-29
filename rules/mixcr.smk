def get_align_target(wildcards):
  if (config["genomic_source"]):
    return "-OvParameters.geneFeatureToAlign=VGeneWithP "
  return ""

rule mixcr_align:
  input:
    os.path.abspath("data/temp/fastq/{sample}_1.fastq"),
    os.path.abspath("data/temp/fastq/{sample}_2.fastq")
  #conda: "../envs/igFinder.yaml"
  log:
    os.path.abspath("logs/mixcr_align/{sample}.log")
  resources:
    mem_mb=16000
  output:
    os.path.abspath("data/mixcr/aligned/{sample}.vdjca")
  threads: 24
  params:
    tgt=get_align_target
  # group: "align"
  shell: 
    "mixcr align -t {threads} -p kAligner2 -f -s hsa "
    "-OallowPartialAlignments=true -OsaveOriginalReads=true "
    "{params.tgt} "
    "{input} {output} 1>{log} 2>&1"


rule mixcr_assemble:
  input:
    rules.mixcr_align.output
  #conda: "../envs/igFinder.yaml"
  output:
    clones=os.path.abspath("data/mixcr/clones/{sample}.clna")
  resources:
    mem_mb=24000
  params:
    Xmx = lambda wildcards, resources: resources.mem_mb - 8000
  log:
    os.path.abspath("logs/mixcr_assemble/{sample}.log")
  # group: "igFinder"
  threads: 8
  shell: 
    "mixcr -Xmx{params.Xmx}m assemble -a -f -t {threads} {input} {output.clones} 1>{log} 2>&1"

rule mixcr_export:
  input:
    rules.mixcr_assemble.output.clones
  #conda: "../envs/igFinder.yaml"
  output:
    short=os.path.abspath("data/mixcr/clone_summary/{sample}_clone_summary.txt"),
    func_clones=os.path.abspath("data/mixcr/func_clones/{sample}_func_clones.txt"),
    full=os.path.abspath("data/mixcr/full_clones/{sample}_full_clones.txt")
  params:
    short = "-cloneId -count -fraction -vGene -dGene -jGene -aaFeature CDR3 -vBestIdentityPercent -vIdentityPercents -jBestIdentityPercent -jIdentityPercents -nFeature CDR3 -avrgFeatureQuality CDR3 -minFeatureQuality CDR3"
  log:
    os.path.abspath("logs/mixcr_export/{sample}.log")
  resources:
    mem_mb=8000
  # group: "igFinder"
  shell:
    "mixcr exportClones -t -o {params.short} "
    "{input} {output.func_clones}; "
    "mixcr exportClones {params.short} "
    "{input} {output.short}; "
    "mixcr exportClones --preset full {input} {output.full}"

rule sig_assemble:
  input:
    clones=rules.mixcr_export.output.func_clones,
    clna=rules.mixcr_assemble.output
  output:
    tempdir = directory(os.path.abspath("tmp/mixcr/{sample}")),
    contigs = os.path.abspath("results/mixcr/assembled/{sample}.vdj.fa")
  resources:
    mem_mb=16000
  log:
    os.path.abspath("logs/mixcr_export_sig_clones/{sample}.log")
  # group: "igFinder"
  threads: 4
  script:
    "../scripts/rev_assemble_clones.py"

rule igblast:
  input:
    contigs = rules.sig_assemble.output.contigs,
    db_V = config["igblast"]["germline_db_V"],
    db_J = config["igblast"]["germline_db_J"],
    db_D = config["igblast"]["germline_db_D"],
    aux_file = config["igblast"]["auxiliary_file"]
  resources:
    mem_mb=8000
  output:
    os.path.abspath("results/igblast/{sample}_igblast_output.txt")
  # group: "igFinder"
  shell:
    "if [ -s {input.contigs} ]\n"
    "then\n"
    "    igblastn -germline_db_V {input.db_V} -germline_db_J {input.db_J} "
    "-germline_db_D {input.db_D} -query {input.contigs} -outfmt 7 "
    "-auxiliary_data {input.aux_file} > {output} \n"
    "else\n"
    "    touch {output}; exit\n"
    "fi" 