library(ggplot2)
library(dplyr)
library(gridExtra)
library(fuzzyjoin)
library(plyr)
library(readr)
library(tidyverse)


cs <- c("#FBA72AFF", "#CB7A5CFF", "#5785C1FF")


# bars <- c("coat", "soap","sanger")
# ana <- "/Users/kalebzuckerman/Desktop/paper_pushing/combined_analysis"
# names <- c("COATswga", "swga2.0","Sanger")
# sz <- "_1GB_"
folder <- "/Users/kalebzuckerman/Desktop/paper_pushing/P2COATswga_Ethiopiasamples/analysis/"
names <- c("Sample 1", "Sample 2", "Sample 3")

f_maps <- read.csv("/Users/kalebzuckerman/Desktop/paper_pushing/P2COATswga_Ethiopiasamples/analysis/maps_sample.csv", header=FALSE)
mp <- data.frame(f_maps) 

fig_2A <- ggplot(mp, aes(x=V1, y=V2, fill=V1)) +
  labs(x="Field Sample", y=expression("Reads Mapped to Pf (%)"), tag="A") +
  scale_fill_manual(values = cs) +
  geom_bar(stat="identity") +
  theme_bw(base_size = 24) +
  theme(legend.position="none",
        axis.text.x = element_text(angle = -30, hjust=0.3, vjust=0.5),
        plot.tag = element_text(size = 36),
        plot.tag.position = c(-0.1, 0.98),
        axis.title.x = element_blank()
  )
fig_2A
# ggsave(filename="fig2A.jpeg", plot=fig_2A, units="px", width=1500, height=1500, dpi="print")


dat <- data.frame(matrix(ncol=3, nrow=0))
file_suff = "covPercs.txt"
for (b in bars) {
  # file_path <- paste(ana, "/barcode", b, "/b", b, sz, file_suff, sep="")
  # file_path <- paste(ana, "/b", b, sz, file_suff, sep="")
  file_path <- Sys.glob(paste0(folder, "/**/*", b, "*_covPercs.txt"))
  f <- read.csv(file_path, header=FALSE)
  f$V3 <- names[which(bars == b)]
  dat <- rbind(dat, f)
}


fig_2C <- ggplot(dat, aes(x=V1, y=V2, color=V3)) + 
  geom_line(linewidth = 2, key_glyph = "point") +
  theme_bw(base_size = 24) + 
  labs(tag="C", color="Sample") +
  scale_color_manual(values = cs) +
  scale_x_continuous(name ="Depth of Coverage", breaks=c(1, 5, seq(10,50, by=10)), expand=c(.01,0), limits=c(0,51)) +
  scale_y_continuous( name = "Percent of Genome Covered", breaks=seq(0,100, by = 20), limits = c(0,100)) +
  geom_vline(xintercept = 10, linetype = "dashed", color="grey") + 
  geom_vline(xintercept = 5, linetype = "dashed", color="grey") +
  guides(color = guide_legend(override.aes = list(shape = 15, size = 8, linetype = 0))) +
  theme(
    plot.tag=element_text(size=36),
    legend.key.spacing.y = unit(0.75, "lines"),
    legend.text = element_text(size = rel(1.3)),
    legend.title = element_text(size = rel(1.6)))
fig_2C

folder_path <- "/Users/kalebzuckerman/Desktop/paper_pushing/P2COATswga_Ethiopiasamples/"

# 1. Glob the folder for files starting with "bmerge"
files <- list.files(path = folder_path, pattern = "^bmerge", full.names = FALSE)

# 2. Create the dataframe and extract sizes
file_df <- data.frame(filename = files) %>%
  extract(
    filename, 
    into = c("sample_label", "size_val", "unit"), 
    regex = "^bmerge_(.*?)_([0-9.]+)(GB|MB)",
    remove = TRUE
  ) %>%
  mutate(
    size_val = as.numeric(size_val),
    size_mb = if_else(unit == "GB", size_val * 1000, size_val)
  ) %>%
  select(sample_label, size_mb) %>%
  filter(sample_label %in% good) %>% 
  mutate(sample_label = names[which(bars == sample_label)])

# View the clean table
print(head(file_df))


fig_2B <- ggplot(file_df, aes(x=sample_label, y=size_mb, fill=sample_label)) +
  labs(x="Field Sample", y=expression("Number of Bases Sequenced (Mb)"), tag="B") +
  scale_fill_manual(values = cs) +
  # scale_y_continuous(breaks=c(0, 500, 1000, 1500)) +
  coord_cartesian(ylim = c(0, 1470)) +
  geom_bar(stat="identity") +
  theme_bw(base_size = 24) +
  theme(legend.position="none",
        axis.text.x = element_text(angle = -30, hjust=0.3, vjust=0.5),
        plot.tag = element_text(size = 36),
        plot.tag.position = c(-0.1, 0.98),
        axis.title.x = element_blank()
  )
fig_2B


# ggsave(file=paste("/Users/kalebzuckerman/Desktop/paper_pushing/final_figs","/fig2_fieldSamples_final.jpeg", sep=""), plot = g, units="px", width=6000, height=4000, dpi="print")


name_lookup <- list(GAAB092413_675V = "Sample 1", GM331_580Y = "Sample 2", SND1071814_441L = "Sample 3")
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
  theme_minimal(base_size = 24) +
  labs(x = "Genomic Position",
       y = "Depth of Coverage",
       tag = "D") +
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
  scale_y_continuous(breaks = c(0, 20, 40)) +
  theme(strip.text = element_text(face = "bold"),
        legend.position = "none",
        # axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.tag = element_text(size = 36),
        plot.margin = margin(t= 5.5, b=15, r=5.5, l=5.5),
        panel.spacing.x = unit(2, "lines"),
        panel.spacing.y = unit(1, "lines")
        # axis.title = element_text(size = 12)
  )

boxes <- ggplot(coverage_data, aes(x = sample_display_name, y = depth, fill = sample_display_name, weight=end-start)) +
  geom_boxplot(outlier.shape = NA) + # Boxplots with transparency
  scale_fill_manual(values = cs) +
  labs(y = "Depth of Coverage", fill = "Sample", tag="E") +
  facet_wrap(~gene, scales = "free", nrow = 1) +
  theme_minimal(base_size = 24) + 
  theme(axis.title.x = element_blank(),
        # axis.text = element_text(size = 12),
        strip.text = element_text(face = "bold"),
        plot.tag = element_text(size = 36),
        plot.margin = margin(t=0, b=30, r=5.5, l=5.5),
        axis.text.x = element_blank(),
        legend.key.spacing.y = unit(0.5, "lines")
        # axis.title.y = element_text(size = 12)
  )

g <- grid.arrange(
  arrangeGrob(
    fig_2A + theme(plot.margin = margin(r = 50, l = 50, t=15, b=15)), 
    fig_2B + theme(plot.margin = margin(r = 50, l = 50, t=15, b=15)), 
    nrow=2, ncol=1
  ), 
  fig_2C, 
  ncol=2, widths=c(2,5))

sav <- grid.arrange(ridges, boxes, nrow=3, heights = c(5,2, 0.25))

final <- grid.arrange(g, sav, nrow=2)


ggsave("/Users/kalebzuckerman/Desktop/paper_pushing/final_figs/fig6_fieldSamples_final.jpg", plot = final, height = 6500, width=6000, units = "px")


