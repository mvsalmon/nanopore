require(tidyverse)

args = commandArgs(trailingOnly = TRUE)
#1 bed file with per-base depth for each target generated from bedtools coverage
#2 run name

# Coverage Function -------------------------------------------------------

coverage <- function(cov_file, run_name, save_dir){
  #calculate depth & coverage for each target in coverage file
  
  cov_data <- cov_file %>%
    group_by(name)
  
  #count bases with non-zero coverage
  coverage <- cov_data %>%
    filter(depthAtPos != 0) %>%
    tally(name = 'coveredBases') %>%
    separate(name, into = c('Transcript_ID', 'HGNC_gene'), sep = '_')
  
  #depth and coverage calculations
  depth_table <- cov_data %>%
    summarise(featureLength = n(),
              meanDepth = mean(depthAtPos),
              ) %>%
    #faster to separate summarised data twice than input once
    separate(name, into = c('Transcript_ID', 'HGNC_gene'), sep = '_') %>%
    #join coverage data
    left_join(coverage) %>%
    #calculate fractional coverage of target
    mutate(fractionalCoverage = coveredBases/featureLength)

  #plot depth per target
  depth_plot <- ggplot(depth_table, aes(x = HGNC_gene, y = meanDepth, label = round(meanDepth, 2))) + 
    geom_bar(stat = 'identity') +
    geom_label() +
    labs(title = sprintf("%s mean depth", run_name),
         y = 'Mean Depth')
  
  ggsave(sprintf("%s_depth_plot.pdf", run_name), depth_plot)
  
  #plot fractional coverage per target
  coverage_plot <- ggplot(depth_table, aes(x = HGNC_gene, y = fractionalCoverage, label = round(fractionalCoverage, 2))) +
    geom_bar(stat = 'identity') +
    geom_label() +
    labs(title = sprintf("%s target coverage", run_name),
         y = 'Coverage')
  
  ggsave(sprintf("%s_coverage_plot.pdf", run_name), coverage_plot)
  }


for(d in dir()){
  cov_file <- read.table(args[1], col.names = c("chrom", "chromStart", "chromEnd", "name", "basePos", "depthAtPos"))
  
  coverage(cov_file, args[2], d)
  }


quit(save = "no")







