require(tidyverse)

args = commandArgs(trailingOnly = TRUE)
#1) run name

# Coverage Function -------------------------------------------------------

coverage <- function(cov_file, run_name, save_dir){
  #calculate depth & coverage for each target in coverage file

  cov_data <- cov_file %>%
    group_by(gene)

  #count bases with non-zero coverage
  coverage <- cov_file %>%
    filter(depthAtPos != 0) %>%
    tally(name = 'coveredBases')
    
    # if(str_detect(coverage$gene, "^ENST")){
    #   coverage <- coverage %>%
    #     separate(gene, into = c('Transcript_ID', 'gene'), sep = '_')
    # }

  #depth and coverage calculations
  depth_table <- cov_data %>%
    summarise(featureLength = n(),
              #mean of per-base depth for a feature
              meanDepth = mean(depthAtPos),
              ) %>%
    #faster to separate summarised data twice than input once
    #separate(gene, into = c('Transcript_ID', 'gene'), sep = '_') %>%
    #join coverage data
    left_join(coverage) %>%
    #calculate fractional coverage of target. 1 = all bases in target present in at least 1 read
    mutate(fractionalCoverage = coveredBases/featureLength)

  #plot depth per target
  depth_plot <- ggplot(depth_table, aes(x = gene, y = meanDepth, label = round(meanDepth, 2))) +
    geom_bar(stat = 'identity') +
    geom_label() +
    labs(title = sprintf("%s mean depth", run_name),
         y = 'Mean Depth') +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

  ggsave(sprintf("%s_depth_plot.pdf", run_name), depth_plot, path = save_dir)

  #plot fractional coverage per target
  coverage_plot <- ggplot(depth_table, aes(x = gene, y = fractionalCoverage, label = round(fractionalCoverage, 2))) +
    geom_bar(stat = 'identity') +
    geom_label() +
    labs(title = sprintf("%s target fractional coverage", run_name),
         y = 'Coverage') +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

  ggsave(sprintf("%s_coverage_plot.pdf", run_name), coverage_plot,  path = save_dir)
  }


for(d in dir()){
  #cov_file <- read.table(args[1], col.names = c("chrom", "chromStart", "chromEnd", "name", "basePos", "depthAtPos"))
  cov_file <- read.table(list.files(d, pattern = "\\.tsv$", full.names = TRUE),
                         col.names = c("chrom", "chromStart", "chromEnd", "gene", "basePos", "depthAtPos")) %>%
    group_by(gene)

  coverage(cov_file, args[2], d)
  }


quit(save = "no")
