#!/bin/bash
samples_path='/Users/gday/Desktop/All_Custom_Capture_BAM'
samples2_path="/farmshare/user_data/gday/mayo/data/BAM_files_Yan"
samples3_path="/farmshare/user_data/gday/mayo/data/Position_Sorted_BAM_BWA"
out_file="20180626_samples.tsv"
echo -e 'sample\tbam' > $out_file
for bamfile in $( ls -P $samples_path/*.bam); do
  basename="${bamfile##*/}"
  basename="${basename%%.bam}"
  # basename="${basename%%_position_sorted*}"
  # basename="${basename%%.sorted.bam}"
  echo -e "${basename}\t${bamfile}" >> $out_file
done
samples2_path="/farmshare/user_data/gday/mayo/data/BAM_files_Yan"
