library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library(fuzzyjoin) 
library(gridExtra)
library(gt)

cs <- c("COATswga" = "#EE6100FF", "swga2.0" = "#1BB6AFFF", "Sanger" = "#9093A2FF")

read_bed <- function(file, sample_name) {
  data <- read_tsv(file, col_names = c("chromosome", "start", "depth"))
  data$sample <- sample_name
  return(data)
}

ana <- "/Users/kalebzuckerman/Desktop/paper_pushing/combined_analysis/longs/"
file_suff <- "_1GB_varCov_long.bed"

# Read in the BED files for four different samples
sample1 <- read_bed(paste(ana, "bcoat", file_suff, sep=""), "COATswga")
sample2 <- read_bed(paste(ana, "bsanger", file_suff, sep=""), "Sanger")
sample3 <- read_bed(paste(ana, "bsoap", file_suff, sep=""), "swga2.0")


coverage_data <- bind_rows(sample1, sample2, sample3) %>% mutate(end = start + 1)

genes <- read_tsv("/Users/kalebzuckerman/Desktop/paper_pushing/3D7_vars.bed", col_names = c("chromosome", "gene_start", "gene_end", "gene")) %>%
  select(!contains("X")) %>%
  mutate(lab = paste0(chromosome, ": ", gene))


# Perform interval-based join to associate each coverage interval with the corresponding gene
coverage_data <- genome_inner_join(
  coverage_data, genes,
  by = c("chromosome" = "chromosome", "start" = "gene_start", "end" = "gene_end")
) %>%
  rename(chromosome = chromosome.x) %>%
  select(!chromosome.y) %>%
  mutate(lab = factor(lab, levels = unlist(genes$lab))) %>% 
  arrange(gene)


b <- ggplot(coverage_data, aes(x = sample, y = depth, fill = sample, weight = depth)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.6) + # Boxplots with transparency
  # geom_blank(aes(x=start)) +
  facet_wrap(~lab, scales = "free", nrow = 6) +
  scale_y_continuous(
    n.breaks = 3,
    labels = function(x) {
      round(x)
    }) +
  scale_fill_manual(values = cs) + 
  theme_minimal() + # Modern theme
  theme(axis.text.x = element_blank(),
        axis.title = element_blank(),
        legend.position = "top",
        panel.grid.major.x = element_blank()) # Adjust theme

ggsave("/Users/kalebzuckerman/Desktop/paper_pushing/final_figs/supplementary/fig_var_supp.jpeg", b,
       width = 12, height = 10, units = "in")


vars <- read_tsv("/Users/kalebzuckerman/Desktop/paper_pushing/3D7_vars.bed", col_names = c("Chromosome", "Start", "End", "Name", "X", "Strand")) %>%
  select(!contains("X"))

vars1 <- vars[1:30, ]
vars2 <- vars[31:nrow(vars), ]


var_table <- vars2 |>
  gt() |>
  cols_align("center") |>
  cols_width(everything() ~ px(150)) |>
  tab_style(style = cell_text(weight = "bold"), locations = cells_column_labels()) |>
  tab_style(
    style = cell_borders(sides = "right", color = "#d3d3d3", weight = px(1)),
    locations = cells_body(columns = c("Chromosome", "Start", "End", "Name"))
  ) |>
  fmt_integer(c("Start", "End"))
var_table

  
gtsave(var_table, "/Users/kalebzuckerman/Desktop/paper_pushing/final_figs/supplementary/var_table2.png")
