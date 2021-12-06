require(tidyverse)

#load arguments
#1 adapteive sampling summary fastqfile
#2 run name

args = commandArgs(trailingOnly = TRUE)
#print(args)

adaptive_seq <- read.csv(args[1])
#print(head(adaptive_seq))

#tally adaptive seq descisions
counts <- count(adaptive_seq, decision) %>%
          mutate(pct = n/sum(n)*100)

write.table(counts, "adaptive_seq_stats.txt", sep = " ", row.names = F)

#split by decision
decisions <- adaptive_seq %>%
  split(f = as.factor(.$decision))

#save reads with each decision - python script takes whole adaptive output file, but may be useful to have separated files....
# for(dec in names(decisions)){
#   x <- decisions[[dec]]
#   write.table(x, paste(args[2], dec, "read.txt", sep = "_"),
#               quote = FALSE,
#               sep = "\t",
#               row.names = FALSE)
#                       }

#plot read descisions
read_count_plot <- ggplot(counts, aes(y = n, x = decision, label = round(pct, 2))) +
      geom_bar(stat = 'identity') +
      geom_label() +
      labs(title = paste(args[2]),
            x = "Decision",
            y = "Number of reads")

ggsave(paste0(args[2], "_adaptive_stats.pdf"), read_count_plot)

quit(save = "no")
