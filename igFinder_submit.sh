#!/bin/bash

INPUT_BAM_DIR=$BAM
OUTPUT_BAM_DIR="$STAFF"

for FILE in $INPUT_BAM_DIR/s_28242.bam; do
    FILENAME="${FILE##*/}"
    FILENAME="${FILENAME%.*}"
    DIR="$OUTPUT_BAM_DIR/${FILENAME}"
    if [[ -d "$DIR" ]]; then
        echo "$DIR exists"
    else
        echo "making $DIR"
        mkdir $DIR
    fi
    OUT="$DIR/mixcr_out"
    
#     if [[ -d "$OUT" ]]; then
#         echo "${OUT} exists; moving useful files before deleting"
#         if [[ -e "${DIR}/${FILENAME}_mixcr_out_rna_aligned.vdjca" ]]; then
#             echo "alignments already done"
#             mv "${DIR}/${FILENAME}_mixcr_out_rna_aligned.vdjca" "${DIR}"
#         fi
#         rm -r $OUT
#     fi

#     mkdir $OUT
    
#     qsub -v DIR=$DIR,FILE=$FILE,OUT=$OUT,FILENAME=$FILENAME \
#     -N ${FILENAME}_igFinder \
#     -o {OUT}/igFinder_${FILENAME}.output.txt \
#     -wd $OUT \
#     -pe threaded 16 \
#      $STAFF/igFinder/igFinder_master.sh
    
    
    
    # qsub -v DIR=$DIR,IGDATA=$IGDATA,FILE=$FILE,OUT=$OUT,FILENAME=$FILENAME,PERL5LIB=/farmshare/user_data/gday/bin/perl5:$PERL5LIB,VELVETOPTIMISER=/farmshare/user_data/gday/VelvetOptimiser-2.2.5/VelvetOptimiser.pl \
    # -N ${FILENAME}_igFinder \
    # -o {OUT}/igFinder_${FILENAME}.output.txt \
    # -wd $OUT \
    # -pe shm 32 \
    # -l large=1 \
    # /farmshare/user_data/gday/scripts/igFinder_master.sh
    
#        --partition=gpu \
    # --gres=gpu:1 \ 
    
    
    
    # sbatch --export=DIR=$DIR,IGDATA=$IGDATA,FILE=$FILE,OUT=$OUT,FILENAME=$FILENAME,PERL5LIB=/farmshare/user_data/gday/bin/perl5:$PERL5LIB,PATH=/farmshare/user_data/gday/bin:/farmshare/user_data/gday/VelvetOptimiser-2.2.5:$PATH,SAMTOOLS=/farmshare/user_data/gday/bin/samtools,VELVETOPTIMISER=/farmshare/user_data/gday/VelvetOptimiser-2.2.5/VelvetOptimiser.pl \
    # --job-name=${FILENAME}_igFinder \
    # --output=${OUT}/igFinder_${FILENAME}.output.txt \
    # --workdir=$OUT \
    # --partition=gpu \
    # --gres=gpu:1 \
    # --cpus-per-task=32 \
    # --mem-per-cpu=4000 \
    # /farmshare/user_data/gday/scripts/igFinder_master.sh
    sleep 10
done
    