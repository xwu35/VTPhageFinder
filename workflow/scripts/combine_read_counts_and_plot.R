#!/usr/bin/env Rscript

# this script is used to combine read counts from different steps and plot the composition of reads

suppressPackageStartupMessages(library("tidyverse"))

# define input
raw_reads <- snakemake@input$raw_reads
after_trimmomatic <- snakemake@input$after_trimmomatic
after_phix <- snakemake@input$after_phix
after_human <- snakemake@input$after_human
after_host_paired <- snakemake@input$after_host_paired
after_host_singleton <- snakemake@input$after_host_singleton

# raw reads
raw_read_counts <- read.table(raw_reads, sep = "\t", header = T) %>% 
  mutate(sample = sub("_R1.*", "", basename(file)),
         total_raw_reads = num_seqs * 2) %>%
  select(sample, total_raw_reads) 

# after trimmomatic
counts_after_trimmomatic <- read.table(after_trimmomatic, 
                                     sep = "\t", header = T) %>% 
  mutate(sample = sub("_R1.*", "", basename(file)),
         total_reads_after_trimmomatic = num_seqs * 2) %>%
  select(sample, total_reads_after_trimmomatic) 

# after phix
counts_after_phix <- read.table(after_phix, 
                                  sep = "\t", header = T) %>% 
  mutate(sample = sub("_R1.*", "", basename(file)),
         total_reads_after_phix = num_seqs * 2) %>%
  select(sample, total_reads_after_phix) 

# after human
counts_after_human <- read.table(after_human, 
                                   sep = "\t", header = T) %>% 
  mutate(sample = sub("_R1.*", "", basename(file)),
         total_reads_after_human = num_seqs * 2) %>%
  select(sample, total_reads_after_human)

# after host
paired_counts_after_host <- read.table(after_host_paired, 
                                  sep = "\t", header = T) %>%
  mutate(sample = sub("_R1.*", "", basename(file)),
         total_paired_reads_after_host = num_seqs * 2) %>%
  select(sample, total_paired_reads_after_host)
 
singleton_counts_after_host <- read.table(after_host_singleton, 
                                        sep = "\t", header = T) %>%
  mutate(sample = gsub(".unmapped.singletons.fastq.gz", "", basename(file))) %>%
  rename(total_singleton_reads_after_host = "num_seqs") %>%
  select(sample, total_singleton_reads_after_host)

# merge all tables
df_lists = list(raw_read_counts, counts_after_trimmomatic, counts_after_phix,
                counts_after_human, paired_counts_after_host, singleton_counts_after_host)

combined_table <- Reduce(function(x, y) merge(x, y, by = "sample", all = TRUE), df_lists)

# calculate number of removed reads at each step
combined_table <- combined_table %>% 
  mutate(adapters_low_quality = total_raw_reads - total_reads_after_trimmomatic,
         phiX = total_reads_after_trimmomatic - total_reads_after_phix,
         human = total_reads_after_phix - total_reads_after_human,
         host = total_reads_after_human - total_paired_reads_after_host - total_singleton_reads_after_host,
         reads_for_assembly = total_paired_reads_after_host + total_singleton_reads_after_host)

# barplot of removed reads at each step
# define colors
colors <- c("#009E73", "#E69F00", "#56B4E9", "#999999", "#CC79A7")

# plot
p_composition_reads <- combined_table %>%
  select(sample, adapters_low_quality, phiX, human, host, reads_for_assembly) %>%
  pivot_longer(!sample, names_to = "type", values_to = "count") %>% 
  ggplot(aes(x=sample, y=count, fill=type)) +
  geom_bar(position="fill", stat="identity") +
  theme_classic() +
  theme(axis.title = element_text(size = 11),
        axis.line = element_line(colour = "black"), 
        axis.text.x=element_text(size=8, colour = "black", angle=90, vjust=0.5),
        axis.text.y = element_text(size=10, colour = "black"),
        axis.ticks.x=element_blank(), 
        legend.title=element_blank(),
        legend.text=element_text(size=8),
        legend.position = "bottom") + 
  guides(fill=guide_legend(nrow=1)) +
  labs(x="", y="Proportion of reads") +
  scale_y_continuous(expand = c(0,0)) +
  scale_fill_manual(values = colors) # + coord_flip()  

# select few columns to write out
output_table <- combined_table %>%
  select(sample, total_raw_reads, adapters_low_quality, phiX, human, host, reads_for_assembly) 
  
# save table
write.table(output_table, snakemake@output$table,
            sep="\t", row.names = F, col.names=T, quote = F)

# save figure
ggsave(file=snakemake@output$figure,
       plot=p_composition_reads, 
       width=9, height=5)
