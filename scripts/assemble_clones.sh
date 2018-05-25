#!/bin/bash

clns_txt=$(realpath $1)
vdjca_file=$(realpath $2)
mixcr_index=$(realpath $3)
OUT=$(realpath $4)
VDJ_seqs=$(realpath $5)
threads=$6
FILENAME=$7


sig_clones=($(awk  'FNR > 1 {if (($3 >= 0.05 || $2 >= 5))  print $1}' $clns_txt))
# will try to assemble all clones with EITHER:
#    1) more than 5 reads
#    2) Clonal fraction over 0.05


if [[ -e "$VDJ_seqs" ]]; then
    echo "erasing past VDJ_seqs file"
    rm -f "$VDJ_seqs"
    echo "making new VDJ_seqs file"
fi

if  [[ -d "$OUT" ]]; then
    rm -rf $OUT
fi

mkdir $OUT


if ! [[ -d "${OUT}/vOpt_temp" ]]; then
    mkdir "${OUT}/vOpt_temp" 
fi

touch "$VDJ_seqs"

for X in "${sig_clones[@]}"; do
    echo "processing clone $X"
    mkdir ${OUT}/vOpt_temp/${X}
    mixcr exportReadsForClones $mixcr_index $vdjca_file ${X} ${OUT}/reads.fastq
    cd ${OUT}/vOpt_temp/${X}
    
    # $VELVETOPTIMISER -v -s 51 -e 151 -x 4 -t $threads -c n50*tbp \
    VelvetOptimiser.pl -v -s 51 -e 151 -x 4 -t $threads -c n50*tbp \
    -k n50*tbp -a -d ${OUT}/cln_${X}_vOptOut   \
    -f " -fastq -shortPaired -separate ${OUT}/reads_cln${X}_R1.fastq ${OUT}/reads_cln${X}_R2.fastq "
    
    if [[ -e "${OUT}/cln_${X}_vOptOut/contigs.fa" ]]; then
        echo "cln ${X} assembly successful"
        awk -v filename=$FILENAME -v x=$X \
        '{ print ((index($0,">")>0))? ">" filename "_clone_" x "_" substr($0,2) : $0 }' \
         "${OUT}/cln_${X}_vOptOut/contigs.fa" >> "$VDJ_seqs"
    else
        echo "cln ${X} assembly failed"
    fi
done