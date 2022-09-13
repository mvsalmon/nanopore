library(tidyverse)

#import mosdepth output *.regions.bed
read_mosdepth <- function(file_path){
  
 depth_df <- read_table(file_path, 
                        col_names = c("chr", "start", "end", "name", "depth"))
 
 return(depth_df)
  
}
#join
K562_all <- bind_rows(lst(K562_run1, K562_run2, K562_run3, K562_run4), .id = "run")

K562_all$name <- str_replace(K562_all$name, "BCR_ABL1_BREAKPOINT", "ABL1 Upstream")

#save combined data
write.table(K562_all, "K562_all_runs_depth.tsv", sep = "\t", quote = F, row.names = F)

#plot depth of k562 runs

depth_plot <- function(depth_df){

  plt <- ggplot(depth_df, aes(x = name, y = depth, fill = run)) + 
  geom_col(position = "dodge") +
  scale_y_continuous(n.breaks = 10) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                   face = 'italic', size = 10)) +
  labs(x = 'Gene', y = 'Depth')

  return(plt)
}

#summary table of depth
K562_depth_summary_table <- K562_all %>% 
  pivot_wider(names_from = run, values_from = depth) %>%
  rename("Run1" = "K562_run1",
         "Run2" = "K562_run2",
         "Run3" = "K562_run3",
         "Run4" = "K562_run4",
         "Gene" = "name") %>%
  select(Gene:Run4) %>%
  arrange(Gene)

write.table(K562_depth_summary_table, file = "K562_depth_summary.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)


#summary stats
K562_all <- K562_all %>% group_by(run)

K562_summary_stats <- K562_all %>%
  summarise(mean_depth = mean(depth), median_depth = median(depth),
            min_depth = min(depth), max_depth = max(depth))

#get read len frequencies
off_target_len <- read_table("off_target_len.txt", 
                              col_names = c('len', 'freq'))

on_target_len <- read_table("on_target_len.txt", 
                              col_names = c('len', 'freq'))
#join
read_lengths <- bind_rows(lst(on_target_len, off_target_len), .id = "location") %>%
  group_by(location)

#plot

read_lengths %>% filter(len < 2000) %>%
ggplot( aes(x = len, y = freq, fill = location)) + 
  labs(y = "Frequency", x = "Read Length",
       fill = "") +
  scale_fill_discrete(labels = c("Off Target", "On Target")) +
  theme(axis.text = element_text(size = 12, face = "bold"),
        text = element_text(size = 12, face = "bold")) +
  geom_col(alpha = 0.5, 
           position = 'identity') #avoids a stacked plot

#% of on target reads above diff lengths for adaptive sampling
#define read length thresholds
threshold <- c(0, 250, 500, 750, 1000)

#get total reads
tot_reads <- sum(read_lengths$freq)

for(t in threshold){
  t_reads <- read_lengths %>% 
    filter(len > t)
   
 on_reads <- filter(t_reads, location == "on_target_len")
 n_on <- sum(on_reads$freq)

 off_reads <- filter(t_reads, location == "off_target_len")
 n_off <- sum(off_reads$freq)*100
 
  fract_off <- n_off/n_on
  print(sprintf("Percentage of reads of len >= %f off target: %f", t, fract_off))

}

##depth plots

depth_df <- separate(col ="name", into = c('ENST', 'HGNC'), sep = "_")
