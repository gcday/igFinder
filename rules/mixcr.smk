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