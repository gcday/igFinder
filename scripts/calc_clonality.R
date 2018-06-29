library(dplyr)


gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

cc <- read.table(snakemake@input[[1]], sep='\t', header=TRUE, fill=TRUE)


cc[is.na(cc)] <- 0

cc_avg <- select(cc, sample, cloneFraction, cloneCount) %>% 
  group_by(sample) %>% 
  summarize(meanFreq = mean(cloneFraction), 
            geomeanFreq = gm_mean(cloneFraction), 
            ratioMeans = mean(cloneFraction) / gm_mean(cloneFraction),
            read_count = sum(cloneCount))

cc_avg <- mutate(cc_avg, clonal = meanFreq >= 0.2 | ratioMeans > 1.25)

write.table(cc_avg, file = snakemake@output[[1]], 
            quote=FALSE, sep = "\t", row.names = FALSE)
