def get_bam(wildcards):
    return samples.loc[(wildcards.sample), ["bam"]].dropna()

rule samtools_filter:
    input:
        bam=get_bam,
        ig_bed=config["ig_bed"]
    output:
        merged=temp("data/filtered_bam/{sample}.bam"),
        temp1=temp("data/filtered_bam/temp_1_{sample}.bam"),
        temp2=temp("data/filtered_bam/temp_2_{sample}.bam")
    log:
        "logs/samtools_filter/{sample}.log"
    threads: 16
    shell:
        "samtools view -uh -@ {threads} -f 13 -o {output.temp1} {input.bam}  && "
        "samtools view -uh -@ {threads} -f 1 -M -L {input.ig_bed} -o {output.temp2} {input.bam} && "
        "samtools merge {output.merged} {output.temp1} {output.temp2} "

        # "samtools view -uh -@ {threads} -f 13 {input.bam} | samtools sort -@ {threads} -o {output.temp1} && "
        # "samtools view -uh -@ {threads} -f 1 -M -L {input.ig_bed} {input.bam} | samtools sort -@ {threads} -o {output.temp2} && "
        # "samtools merge {output.merged} {output.temp1} {output.temp2} "

rule samtools_sort:
    input:
        "data/filtered_bam/{sample}.bam"
    output:
        temp("data/filtered_namesorted_reads/{sample}.bam")
    log:
        "logs/samtools_sort/{sample}.log"
    threads: 16
    shell:
        "samtools sort -t /tmp -n -m 2000M -@ {threads} -o {output} {input}"

rule samtools_fastq:
    input:
        "data/filtered_namesorted_reads/{sample}.bam"
    output:
        one=temp("data/fastq/{sample}_1.fastq"),
        two=temp("data/fastq/{sample}_2.fastq"),
        singleton=temp("data/fastq/{sample}_singleton.fastq")
    log:
        "logs/samtools_fastq/{sample}.log"
    threads: 16
    shell:
        "samtools fastq -N -@ {threads} "
        "-1 {output.one} -2 {output.two} -s {output.singleton} {input}"
