#!/bin/bash

# requirements: mixcr must be in your $PATH
#               $VELVETOPTIMISER must point to VelvetOptimiser.pl (tested with 2.2.5)

namesort_bam="$DIR/${FILENAME}.namesorted.bam"

input_fastq_1="$DIR/${FILENAME}_1.fastq"

input_fastq_2="$DIR/${FILENAME}_2.fastq"

input_fastq_singleton="$DIR/${FILENAME}_singleton.fastq"

vdjca_file="${DIR}/${FILENAME}_mixcr_out_rna_aligned.vdjca"

VDJ_seqs="${OUT}/${FILENAME}_VDJ_seqs.fa"

cores=`nproc` # default behavior: use all available cores for
              # multithreaded applications
echo $namesort_bam

mixcr

if [[ -e "$vdjca_file" ]]; then
    echo "moving straight to mixcr assembly step"
        
else
    if [[ -e "$input_fastq_1" && -e "$input_fastq_2" ]]; then
        echo "fastq files already present"
    else
        if [[ -e "$namesort_bam" ]]; then
            echo "namesorted bam exists"
        else
            echo "sorting bam"
            samtools sort -@ $cores -n -o "$namesort_bam" $FILE
        fi
        echo "making fastq files from bam"
        # UNLESS YOU ABSOLUTELY KNOW THE INPUT BAM FILE IS NAMESORTED, 
        # DO NOT SKIP THIS STEP!!!
        # SAMTOOLS FASTQ OUTPUTS IN THE SAME ORDER AS IN THE ORIGINAL BAM
        samtools fastq \
        -N \
        -@ $cores \
        -1 "$input_fastq_1" \
        -2 "$input_fastq_2" \
        -s "$input_fastq_singleton" \
        "$namesort_bam"
    fi
    
    
    
    
    
    mixcr align -t $cores -g -a -f -p rna-seq -s hsa \
    -OvParameters.geneFeatureToAlign=VGeneWithP \
    -OallowPartialAlignments=true \
    "$input_fastq_1" \
    "$input_fastq_2" \
    "$vdjca_file"
    
    # -t $cores : Sets number of threads for running in parallel 
    # -p rna-seq : This script is designed to process genomic DNA reads, not RNA-seq.
    #             However, MiXCR reccomends using the RNA-seq settings in order to better
    #             rescue reads covering only portions of a rearranged IG gene.
    # -f : forces overwriting of output files if they exist
    # -g, -a : saves original read information, needed for exporting
    # -s hsa : species=human 
    # -OvParameters.geneFeatureToAlign=VGeneWithP : aligns to germline DNA sequence (not spliced sequence)
    # -OallowPartialAlignments=true : partial alignments help when a read covers either the V segment and
    #                                 part of the CDR3 or the J segment plus part of the CDR3 but not both
fi


# mixcr assemblePartial ${OUT}/${FILENAME}_mixcr_out_rna.3.vdjca ${OUT}/${FILENAME}_alignments_rescued_1.vdjca
# mixcr assemblePartial ${OUT}/${FILENAME}_alignments_rescued_1.vdjca ${OUT}/${FILENAME}_alignments_rescued_2.vdjca
# mixcr assemble -f ${OUT}/${FILENAME}_alignments_rescued_2.vdjca ${OUT}/${FILENAME}_alignments_rescued_clones.clns


mixcr_index="${OUT}/${FILENAME}_reads_by_clone.index"
clns_file="${OUT}/${FILENAME}_clones.clns"
clns_txt="${clns_file}.txt"

mixcr assemble -f -i $mixcr_index $vdjca_file $clns_file

mixcr exportClones -cloneId -count -fraction -vGene -dGene -jGene \
-aaFeature CDR3 -vBestIdentityPercent -vIdentityPercents  \
-jBestIdentityPercent -jIdentityPercents -nFeature CDR3 \
-avrgFeatureQuality CDR3 -minFeatureQuality CDR3 \
$clns_file $clns_txt

# mkdir ${OUT}/vOpt_temp

sig_clones=($(awk  'FNR > 1 {if (($3 >= 0.05 || $2 >= 5))  print $1}' $clns_txt))
# will try to assemble all clones with EITHER:
#    1) more than 5 reads
#    2) Clonal fraction over 0.05

if [[ -e "$VDJ_seqs" ]]; then
    echo "erasing past VDJ_seqs file"
    rm -f "$VDJ_seqs"
else
    echo "making new VDJ_seqs file"
fi

touch "$VDJ_seqs"

for X in "${sig_clones[@]}"; do
    echo "processing clone $X"
    mkdir ${OUT}/vOpt_temp/${X}
    $MIXCR exportReadsForClones $mixcr_index vdjca_file ${X} ${OUT}/reads.fastq
    cd ${OUT}/vOpt_temp/${X}
    
    $VELVETOPTIMISER -v -s 51 -e 151 -x 4 -t 8 -c n50*tbp \
    -k n50*tbp -a -d ${OUT}/cln_${X}_vOptOut_${FILENAME} \
    -f " -fastq -shortPaired -separate ${OUT}/reads_cln${X}_R1.fastq ${OUT}/reads_cln${X}_R2.fastq "
    
    
    if [[ -e "${OUT}/cln_${X}_vOptOut_${FILENAME}/contigs.fa" ]]; then
        echo "cln ${X} assembly successful"
        awk -v filename=$FILENAME -v x=$X \
        '{ print ((index($0,">")>0))? ">" filename "_clone_" x "_" substr($0,2) : $0 }' \
        ${OUT}/cln_${X}_vOptOut_${FILENAME}/contigs.fa >> "${OUT}/${FILENAME}_VDJ_seqs.fa"
    else
        echo "cln ${X} assembly failed"
    fi
done

    


# mixcr exportReadsForClones ${OUT}/${FILENAME}_reads_by_clone.index ${OUT}/${FILENAME}_mixcr_out_rna.3.vdjca 0 ${OUT}/reads.fastq

# srun $VELVETOPTIMISER -v -s 15 -e 151 -x 4 -t 48 -c max -a -d ${OUT}/vOptOut_${FILENAME} -f " -fastq -shortPaired -separate ${OUT}/reads_cln0_R1.fastq ${OUT}/reads_cln0_R2.fastq "

rm -f "$DIR/${FILENAME}_1.fastq"
rm -f "$DIR/${FILENAME}_2.fastq"
rm -f "$DIR/${FILENAME}.namesorted.bam"
rm -f "$DIR/${FILENAME}_singleton.fastq"

exit
