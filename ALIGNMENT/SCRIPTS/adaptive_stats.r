require(tidyverse)

#load arguments
#1 adaptive sampling summary fastqfile
#2 run name


args = commandArgs(trailingOnly = TRUE)


# Adaptive Stats function -------------------------------------------------


adaptive_stats <- function(adaptive_seq, run_name){
#descriptive stats for nanopore adapitive sampling summary file
  
  #tally adaptive seq descisions
  counts <- count(adaptive_seq, decision) %>%
            mutate(pct = n/sum(n)*100)

  write.table(counts, sprintf("%s_adaptive_seq_stats.txt", run_name), sep = " ", row.names = F, quote = F)

  #split by decision
  decisions <- adaptive_seq %>%
    split(f = as.factor(.$decision))


  #plot read descisions
  read_count_plot <- ggplot(counts, aes(y = n, x = decision, label = round(pct, 2))) +
        geom_bar(stat = 'identity') +
        geom_label() +
        labs(title = paste(run_name),
             subtitle = 'Boxes show % of reads in each group',
              x = "Decision",
              y = "Number of reads") +
    theme(axis.text.x = element_text(angle = 45))

  ggsave(sprintf("%s_adaptive_stats.pdf", run_name), read_count_plot)
  }

adaptive_stats(read.csv(args[1]), args[2])

quit(save = "no")
