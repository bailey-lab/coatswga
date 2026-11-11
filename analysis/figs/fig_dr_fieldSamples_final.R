library(ggplot2)
library(dplyr)
library(gridExtra)
library(fuzzyjoin)
library(plyr)
library(readr)

cs <- c("#FBA72AFF", "#CB7A5CFF", "#5785C1FF")

# names <- c("Sample 1", "Sample 2", "Sample 3")

file_paths <- Sys.glob("/Users/kalebzuckerman/Desktop/paper_pushing/P2COATswga_Ethiopiasamples/analysis/*/*_DRCov.bed")


coverage_data <- map_df(file_paths, function(f) {
  id_from_file <- str_extract(basename(f), "(?<=bmerge_)(.*?)(?=_[0-9]+\\.)")
  read_tsv(f, col_names = c("chromosome", "start", "end", "depth"), show_col_types = FALSE) %>%
    mutate(
      sample = id_from_file,
      sample_display_name = as.character(ifelse(id_from_file %in% names(name_lookup), 
                                   name_lookup[id_from_file], 
                                   id_from_file))
    )
})

# Convert chromosome column to factor for proper ordering
coverage_data$chromosome <- factor(coverage_data$chromosome, levels = unique(coverage_data$chromosome))

# Read gene annotation file
genes <- read_tsv("/Users/kalebzuckerman/Desktop/paper_pushing/3D7_DR_genes_RefSeq.bed", col_names = c("chromosome", "gene_start", "gene_end", "gene"))

coverage_data <- coverage_data %>%
  filter(sample %in% good) %>%
  left_join(genes %>% select(chromosome, gene)) 

# Plot using ggplot2
ridges <- ggplot(coverage_data, aes(xmin = start, xmax = end, ymin = 0, ymax = depth, fill = sample_display_name)) +
  geom_rect(size = 0.1) +
  # facet_grid(sample ~ chromosome, scales = "free_x", space = "free_x") +
  facet_grid(sample_display_name ~ gene, scales = "free", space = "free_x") +
  theme_minimal(base_size = 20) +
  labs(x = "Genomic Position",
       y = "Depth of Coverage",
       tag = "A") +
  coord_cartesian(ylim = c(0, 50)) +
  scale_fill_manual(values = cs) +
  scale_x_continuous(expand = c(0.01,0.01),
                     breaks = function(x) {
                       ind <- which.min(abs(unlist(genes["gene_start"]) - x[1]))
                       start <- round_any(genes[[ind, "gene_start"]]*1.0001, 10, f = ceiling)
                       end <- round_any(genes[[ind, "gene_end"]]*0.9999, 10, f = floor)
                       seq(start, end, by = (end - start) / 2)
                     },
                     labels = function(x) {
                       # Format the labels as the actual genomic positions
                       as.character(round(x))
                     }
  ) +
  theme(strip.text = element_text(face = "bold"),
        legend.position = "none",
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.margin = margin(t= 5.5, b=15, r=5.5, l=5.5),
        panel.spacing = unit(0.5, "lines")
        # axis.title = element_text(size = 12)
  )

boxes <- ggplot(coverage_data, aes(x = sample_display_name, y = depth, fill = sample_display_name, weight=end-start)) +
  geom_boxplot(outlier.shape = NA) + # Boxplots with transparency
  scale_fill_manual(values = cs) +
  labs(y = "Depth of Coverage", fill = "Field Sample", tag="B") +
  facet_wrap(~gene, scales = "free", nrow = 1) +
  theme_minimal(base_size = 16) + 
  theme(axis.title.x = element_blank(),
        axis.text = element_text(size = 12),
        strip.text = element_text(face = "bold"),
        plot.margin = margin(t=0, b=30, r=5.5, l=5.5),
        axis.text.x = element_blank(),
        # axis.title.y = element_text(size = 12)
  )

sav <- grid.arrange(ridges, boxes, nrow=3, heights = c(5,2, 0.25))

ggsave("/Users/kalebzuckerman/Desktop/paper_pushing/final_figs/figDR_field_samples_final.jpg", plot = sav, height = 2500, width=4000, units = "px")
# ggsave(paste0("/Users/kalebzuckerman/Desktop/placement_analysis/", "fig_DR.jpg"), plot = ridges,  width = 10, height = 7, )

