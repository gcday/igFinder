def get_bam(wildcards):
    return samples.loc[(wildcards.sample), ["bam"]].dropna()

rule samtools_filter:
    input:
        bam=get_bam,
        ig_bed=config["ig_bed"]
    conda: "../envs/igFinder.yaml"
    output:
        merged=temp("data/filtered_bam/{sample}.bam"),
        temp1=temp("data/filtered_bam/temp_1_{sample}.bam"),
        temp2=temp("data/filtered_bam/temp_2_{sample}.bam")
    log:
        "logs/samtools_filter/{sample}.log"
    threads: 4
    #group: "igFinder"
    shell:
      "samtools view -uh -@ {threads} -f 13 -o {output.temp1} {input.bam} && "
      "samtools view -uh -@ {threads} -f 1 -M -L {input.ig_bed} " 
      " -o {output.temp2} {input.bam} && "
      "samtools merge {output.merged} {output.temp1} {output.temp2} "
      #"samtools index -@ {threads} {input.bam} && "
rule samtools_sort:
    input:
        "data/filtered_bam/{sample}.bam"
    conda: "../envs/igFinder.yaml"
    output:
        temp("data/filtered_namesorted_reads/{sample}.bam")
    log:
        "logs/samtools_sort/{sample}.log"
    threads: 4
    #group: "igFinder"
    shell:
        "samtools sort -t /tmp -n -m 2000M -@ {threads} -o {output} {input}"

rule samtools_fastq:
    input:
        "data/filtered_namesorted_reads/{sample}.bam"
    conda: "../envs/igFinder.yaml"
    output:
        one=temp("data/fastq/{sample}_1.fastq"),
        two=temp("data/fastq/{sample}_2.fastq"),
        singleton=temp("data/fastq/{sample}_singleton.fastq")
    log:
        "logs/samtools_fastq/{sample}.log"
    threads: 4
    #group: "igFinder"
    shell:
        "samtools fastq -N -@ {threads} "
        "-1 {output.one} -2 {output.two} -s {output.singleton} {input} && "
        "touch {output.singleton}"
rule count_reads:
    input:
        get_bam
    output:
        "data/read_counts/{sample}.txt"
    threads: 4
    shell:
        "samtools view -@ 16 -c {input} > {output}"