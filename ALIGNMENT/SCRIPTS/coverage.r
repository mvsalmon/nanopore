require(tidyverse)

args = commandArgs(trailingOnly = TRUE)
#1 bed file with per-base depth for each target generated from bedtools coverage
#2 run name

# Coverage Function -------------------------------------------------------

coverage <- function(cov_file, run_name){
  #calculate depth & coverage for each target in coverage file
  
  cov_data <- cov_file %>%
    group_by(name)
  
  #depth and coverage calculations
  depth_table <- cov_data %>%
    summarise(featureLength = n(),
              meanDepth = mean(depthAtPos),
              ) %>%
    separate(name, into = c('Transcript_ID', 'HGNC_gene'), sep = '_')
  
  #plot depth per target
  ggplot(depth_table, aes(x = HGNC_gene, y = meanDepth)) + 
    geom_col(stat = 'identity') +
    labs(title = run_name,
         y = 'Mean Depth')
  }


cov_file <- (read.table(args[1], col.names = c("chrom", "chromStart", "chromEnd", "name", "basePos", "depthAtPos")))

coverage(cov_file, args[2])

quit(save = "no")







