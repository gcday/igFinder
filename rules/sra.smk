def get_srr(wildcards):
    return list(samples[samples["sample"] == wildcards.sample]["SRR"])[0]
    # return list(sra_table[sra_table['Sample_Name'].str.contains(str(wildcards.sample))]['Run'])[0]
rule download:
    params:
        srr=get_srr,
        sra_dir = config["sra_dir"],
        log_path = os.path.abspath("logs/download/{sample}.log")
    output:
        fastq1=temp(os.path.abspath("data/temp/fastq/{sample}_1.fastq")),
        fastq2=temp(os.path.abspath("data/temp/fastq/{sample}_2.fastq")),
        # bam_file="data/bam/{sample}.bam"
    log:
        os.path.abspath("logs/download/{sample}.log")
    resources:
        mem_mb=48000
    # group: "sra"
    threads: 24
    shell:
        "cd {params.sra_dir}; fasterq-dump -f -e {threads} -m 12000MB -p {params.srr} 1>{log} 2>&1\n"
        "if [[ $(stat -c%s {params.srr}_1.fastq) -eq $(stat -c%s {params.srr}_2.fastq) ]]; then \n"
        "   mv {params.srr}_1.fastq {output.fastq1} && mv {params.srr}_2.fastq {output.fastq2} \n"
        "else \n"
        "   echo \"Error in fasterq-dump download, please try again.\" >> {log}\n"
        "   exit 1\n"
        "fi\n"
