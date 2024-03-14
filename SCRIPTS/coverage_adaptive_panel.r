library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

#1) bedtools .tsv coverage output file
#2) run name
#3) output directory

# Coverage Function -------------------------------------------------------

coverage <- function(cov_file, run_name, output_dir) {
  #calculate depth & coverage for each target in coverage file

  #summary of overall sequencing depth
  depth_summary <- cov_file %>%
    ungroup() %>%
    summarise(MinDepth = min(depthAtPos),
              MaxDepth = max(depthAtPos),
              MedianDepth = median(depthAtPos),
              MeanDepth = mean(depthAtPos),
              SD = sd(depthAtPos))
  write_delim(depth_summary,
              sprintf("%s/%s_depth_summary.tsv", output_dir, run_name),
              delim = "\t",
              quote = "none")

  #per-gene summary
  per_gene_depth_summary <- cov_file %>%
    group_by(gene, chrom) %>%
    summarise(MinDepth = min(depthAtPos),
              MaxDepth = max(depthAtPos),
              MedianDepth = median(depthAtPos),
              MeanDepth = mean(depthAtPos),
              SD = sd(depthAtPos))
  write_delim(per_gene_depth_summary,
              sprintf("%s/%s_per_gene_depth_summary.tsv", output_dir, run_name),
              delim = "\t",
              quote = "none")

  #count bases with non-zero coverage
  coverage_data <- cov_file %>%
    filter(depthAtPos != 0) %>%
    tally(name = "coveredBases")


  #depth and coverage calculations
  depth_table <- cov_file %>%
    summarise(
      featureLength = n(),
      medianDepth = median(depthAtPos), #median of per-base depth for a feature
    ) %>%
    #join coverage data
    left_join(coverage_data) %>%
    #calculate fractional coverage of target. 
    # 1 = all bases in target present in at least 1 read
    mutate(fractionalCoverage = coveredBases/featureLength) %>%
    #separate gene name and RefSeq accession
    separate(gene, into = c('gene', 'RefSeq'), sep = ',') %>%
    #output table with per-target depth and coverage info
    write_delim(sprintf("%s/%s_coverage_data.tsv", output_dir, run_name),
                delim = "\t",
                quote = "none")

  #summary of coverage data

  #plot depth per target
  depth_plot <- ggplot(depth_table, aes(x = gene, y = medianDepth)) +
    geom_bar(stat = "identity") +
    labs(title = sprintf("%s median depth", run_name),
         y = "Median Depth") +
    theme(axis.text.x = element_text(size = 4)) +
    scale_x_discrete(guide = guide_axis(n.dodge = 3, angle = 45))

  ggsave(sprintf("%s_depth_plot.pdf", run_name), depth_plot, path = output_dir)

  #plot fractional coverage per target
  coverage_plot <- ggplot(depth_table, aes(x = gene, y = fractionalCoverage)) +
    geom_bar(stat = "identity") +
    theme(axis.text.x = element_text(angle = 45, 
                                     hjust = 1, vjust = 0.5, size = 4)) +
    labs(title = sprintf("%s target fractional coverage", run_name),
         y = "Coverage")

  ggsave(sprintf("%s_coverage_plot.pdf", run_name),
         coverage_plot,  path = output_dir)
}

cov_file <- read_delim(
  args[[1]],
  delim = "\t",
  col_names = c("chrom", "chromStart", "chromEnd", "gene",
                "basePos", "depthAtPos")
) %>%
  #grouping by chrom accounts for PAR1 on X/Y where no mapping to Y makes a 
  #few genes in this region appear to only have half the depth/coverage
  #chrom group gives a row for each gene on X and Y - may be useful to keep
  #as internal control?
  group_by(gene)

coverage(cov_file, args[[2]], args[[3]])

quit(save = "no")
