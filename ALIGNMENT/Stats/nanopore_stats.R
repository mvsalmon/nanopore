library(tidyverse)

#import depth bed files
read_mosdepth <- function(file_path){
  
 depth_df <- read_table(file_path, 
                        col_names = c("chr", "start", "end", "name", "depth"))
 
 return(depth_df)
  
}
#join
K562_all <- bind_rows(lst(K562_run1, K562_run2, K562_run3), .id = "run")
#plot depth of k562 runs

depth_plot <- function(depth_df, run_name = NULL, single_run = TRUE){
  if(single_run){
   fill_col <- NULL}
  else{
    fill_col <- "run"
  }
  
  plt <- ggplot(depth_df, aes_(x = quote(name), y = quote(depth), fill = quote(fill_col))) + 
  geom_col(position = "dodge") +
  scale_y_continuous(n.breaks = 10) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 1)) +
  labs(title = sprintf("%s sequencing depth per target", run_name),
       x = 'Gene', y = 'Depth')

  return(plt)
}

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
  labs(title = "K562 library 1 - adaptive sampling") +
  theme(axis.text = element_text(size = 12, face = "bold")) +
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
