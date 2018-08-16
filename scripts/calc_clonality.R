library(dplyr)

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
# setwd("/Users/grady/OneDrive - Leland Stanford Junior University/Mayo")

cc <- read.table(snakemake@input[[1]], sep='\t', header=TRUE, fill=TRUE)
reads_table <- read.table(snakemake@input[[2]], sep='\t', header=TRUE, fill=TRUE)
# cc <- read.table("20180701_custom_capture_clones.tsv", sep='\t', header=TRUE, fill=TRUE)

cc[is.na(cc)] <- 0
cc_with_chains <- mutate(cc, chain = substr(bestVGene, 1, 3)) %>% 
  group_by(sample, chain) %>% arrange( desc(sample), desc(cloneCount)) 


by_chain = summarize(cc_with_chains,
                     top_two = if (n() >= 2) first(cloneCount) + nth(cloneCount, 2) else first(cloneCount))
most_common <- group_by(by_chain, sample) %>% summarize(singleCloneReads = sum(top_two))

cc_avg <- select(cc, sample, cloneFraction, cloneCount) %>% 
  group_by(sample) %>% 
  summarize(meanFreq = mean(cloneFraction), 
            geomeanFreq = gm_mean(cloneFraction), 
            ratioMeans = mean(cloneFraction) / gm_mean(cloneFraction),
            vdj_read_count = sum(cloneCount))

cc_avg <- mutate(cc_avg, clonal = meanFreq >= 0.2 | ratioMeans > 1.25)
cc_frac <- left_join(cc_avg, most_common, by = "sample") %>% mutate(singleCloneFraction = singleCloneReads / vdj_read_count)
cc_frac[is.na(cc_frac)] <- 0
cc_final <- mutate(cc_frac, updated_clonality = singleCloneFraction > 0.5)
cc_final <- left_join(cc_final, reads_table, by = "sample")
write.table(cc_final, file = snakemake@output[[1]], 
            quote=FALSE, sep = "\t", row.names = FALSE)