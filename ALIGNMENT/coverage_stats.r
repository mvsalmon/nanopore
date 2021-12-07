cov_table <- read.table(BEDTOOLS_OUT,
                       col_names = c("chr", "start", "end", "name", "n_reads", "covered_bases", "region_len", "fract_covered"))

ggplot(cov_table, aes(x = name, y = fract_covered)) + geom_col() + theme(axis.text.x = element_text(angle = 90))