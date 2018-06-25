#!/bin/bash
samples_path="/farmshare/user_data/gday/mayo/data/BAM_bwa-mem_position_sorted_yan"
samples2_path="/farmshare/user_data/gday/mayo/data/BAM_files_Yan"

out_file="20180624_samples.tsv"
echo -e 'sample\tbam' > $out_file
for bamfile in $( ls $samples_path/*.bam $samples2_path/*.bam); do
  basename="${bamfile##*/}"
  basename="${basename%%_position_sorted*}"
  echo -e "${basename}\t${bamfile}" >> $out_file
done
samples2_path="/farmshare/user_data/gday/mayo/data/BAM_files_Yan"
