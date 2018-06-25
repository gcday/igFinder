#!/bin/bash

# bam_dir=/farmshare/user_data/gday/igFinder

# OUT=/farmshare/user_data/gday/igFinder/cc_results.txt

# clones_dir=$bam_dir/clones

# if [[ ! -d $clones_dir ]]; then
#     mkdir $clones_dir
# fi

# if [[ -e $OUT ]]; then
#     rm -f $OUT
# fi

# touch $OUT
# i=0

# for bam_file in $bam_dir/*.bam; do
#     sample=${bam_file##*/}
#     sample=${sample%.bam}

#     # echo $sample
    
#     clns_file=$bam_dir/$sample/mixcr_out/${sample}_clones.clns
    
#     cnls_txt_file=${clns_file}.txt
    
#     if [[ ! -e $clns_file ]]; then
#         echo "WARNING: $clns_file missing" 
#     else
        
#         cp $clns_file $clones_dir
#         mixcr exportClones --preset full "$clns_file" "${clones_dir}/${sample}_clones.clns.txt"
        
#         if [ "$i" -eq 0 ]; then
#             awk -v bam_file=$bam_file 'FNR == 1 { print "sample\t",$0 }' $cnls_txt_file >> $OUT
#         fi

#         # echo $DIR
#         let i+=1

#         # # echo $FILE
#         awk -v bam_file=$bam_file 'FNR > 1 { print bam_file,"\t",$0}' $cnls_txt_file >> $OUT
#     fi
# done


out_file=/farmshare/user_data/gday/igFinder/clones/vdjtools/metadata.txt

echo -e "#file.name\tsample.id" > $out_file

{
for vdjtools_file in /farmshare/user_data/gday/igFinder/clones/vdjtools.*_clones.clns.txt; do
    sample_name=${vdjtools_file##*/vdjtools.}
    sample_name=${sample_name%_clones.clns.txt}
    echo -e "${vdjtools_file}\t${sample_name}" >> $out_file
done
}










