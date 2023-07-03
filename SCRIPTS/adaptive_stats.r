library(tidyverse)

#load arguments
#1 adaptive sampling summary file
#2 run name
#3 output dir

args <- commandArgs(trailingOnly = TRUE)

sprintf("INPUT FILE: %s", args[1])
sprintf("RUN NAME: %s", args[2])

# Adaptive Stats function -------------------------------------------------


adaptive_stats <- function(adaptive_seq, run_name, output_dir) {
#descriptive stats for nanopore adapitive sampling summary file
#tally adaptive seq descisions
  counts <- count(adaptive_seq, decision) %>%
            mutate(pct = n/sum(n)*100)

  write.table(counts,
  sprintf("%s/%s_adaptive_seq_stats.txt", output_dir, run_name),
  sep = " ", row.names = FALSE, quote = FALSE)

  #split by decision
  decisions <- adaptive_seq %>%
    split(f = as.factor(.$decision))


  #plot read descisions
  read_count_plot <- ggplot(counts,
        aes(y = n, x = decision, label = round(pct, 2))) +
        geom_bar(stat = 'identity') +
        geom_label() +
        labs(title = paste(run_name),
             subtitle = "Boxes show % of reads in each group",
              x = "Decision",
              y = "Number of reads") +
    theme(axis.text.x = element_text(angle = 45))

  ggsave(sprintf("%s/%s_adaptive_stats.pdf", output_dir, run_name),
  read_count_plot)
  }

summary_file <- read.csv(args[1])
adaptive_stats(summary_file, args[2], args[3])

quit(save = "no")
