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
    resources:
        mem_mb=24000
    # group: "sra"
    threads: 8
    shell:
        "cd {params.sra_dir}; fasterq-dump -p {params.srr} 1>{params.log_path} 2>&1 && "
        "mv {params.srr}_1.fastq {params.fq1_path} && mv {params.srr}_2.fastq {params.fq2_path} "