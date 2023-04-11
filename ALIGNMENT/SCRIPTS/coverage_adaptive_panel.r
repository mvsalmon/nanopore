require(tidyverse)

args = commandArgs(trailingOnly = TRUE)
#1) run name

# Coverage Function -------------------------------------------------------

coverage <- function(cov_file, run_name, save_dir){
  #calculate depth & coverage for each target in coverage file

  #  cov_data <- cov_file %>%
  #    group_by(gene)

  #count bases with non-zero coverage
  coverage_data <- cov_file %>%
    filter(depthAtPos != 0) %>%
    tally(name = 'coveredBases')

    # if(str_detect(coverage$gene, "^ENST")){
    #   coverage <- coverage %>%
    #     separate(gene, into = c('Transcript_ID', 'gene'), sep = '_')
    # }

  #depth and coverage calculations
  depth_table <- cov_file %>%
    summarise(featureLength = n(),
              #mean of per-base depth for a feature
              meanDepth = mean(depthAtPos),
              ) %>%
    #join coverage data
    left_join(coverage_data) %>%
    #calculate fractional coverage of target. 1 = all bases in target present in at least 1 read
    mutate(fractionalCoverage = coveredBases/featureLength) %>%
    #separate gene name and RefSeq accession
    separate(gene, into = c('gene', 'RefSeq'), sep = ',') %>%
    #output table with per-target depth and coverage info
    write_delim(sprintf("%s_coverage_data.tsv", run_name),
                delim = "\t",
                quote = "none")

  #summary of coverage data

  #plot depth per target
  depth_plot <- ggplot(depth_table, aes(x = gene, y = meanDepth)) +
    geom_bar(stat = 'identity') +
    labs(title = sprintf("%s mean depth", run_name),
         y = 'Mean Depth') +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

  ggsave(sprintf("%s_depth_plot.pdf", run_name), depth_plot, path = save_dir)

  #plot fractional coverage per target
  coverage_plot <- ggplot(depth_table, aes(x = gene, y = fractionalCoverage)) +
    geom_bar(stat = 'identity') +
    labs(title = sprintf("%s target fractional coverage", run_name),
         y = 'Coverage') +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

  ggsave(sprintf("%s_coverage_plot.pdf", run_name), coverage_plot,  path = save_dir)
  }


#for(d in dir()){
  #cov_file <- read.table(args[1], col.names = c("chrom", "chromStart", "chromEnd", "name", "basePos", "depthAtPos"))
  cov_file <- read_delim(list.files(pattern = "\\.tsv$", full.names = TRUE),
                        delim = "\t",
                        col_names = c("chrom", "chromStart", "chromEnd", "gene", "score", "strand", "basePos", "depthAtPos")
                        ) %>%
                        
    group_by(gene)

  coverage(cov_file, args[1], getwd())



quit(save = "no")
