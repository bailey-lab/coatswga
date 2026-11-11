library(ggplot2)
library(dplyr)
library(gridExtra)
library(fuzzyjoin)
library(plyr)
library(readr)

cs <- c("COATswga" = "#EE6100FF", "swga2.0" = "#1BB6AFFF", "Sanger" = "#9093A2FF")

bars <- c("coat", "soap","sanger")
ana <- "/Users/kalebzuckerman/Desktop/paper_pushing/combined_analysis"
names <- c("COATswga", "swga2.0","Sanger")
sz <- "_1GB_"


f_maps <- read.csv("/Users/kalebzuckerman/Desktop/paper_pushing/final_figs/maps.csv", header=FALSE)
mp <- data.frame(f_maps)

fig_2A <- ggplot(mp, aes(x=V1, y=V2, fill=V1)) +
  labs(x="Primer Set", y=expression("Percent of Reads Mapped to Pf"), tag="A") +
  scale_fill_manual(values = cs) +
  geom_bar(stat="identity") +
  theme_bw(base_size = 24) +
  theme(legend.position="none", 
        plot.tag = element_text(size=36),
        axis.text.x = element_text(angle = 30, hjust=1))
fig_2A
# ggsave(filename="fig2A.jpeg", plot=fig_2A, units="px", width=1500, height=1500, dpi="print")


dat <- data.frame(matrix(ncol=3, nrow=0))
file_suff = "covPercs.txt"
for (b in bars) {
  # file_path <- paste(ana, "/barcode", b, "/b", b, sz, file_suff, sep="")
  file_path <- paste(ana, "/b", b, sz, file_suff, sep="")
  f <- read.csv(file_path, header=FALSE)
  f$V3 <- names[which(bars == b)]
  dat <- rbind(dat, f)
}

fig_2B <- ggplot(dat, aes(x=V1, y=V2, color=V3)) + 
  geom_line(linewidth = 2, key_glyph = "point") +
  theme_bw(base_size = 24) + 
  labs(tag="B", color="Primer Set") +
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
fig_2B

# ggsave(file=paste("/Users/kalebzuckerman/Desktop/paper_pushing/final_figs","/fig2.jpeg", sep=""), plot = g, units="px", width=4400, height=2120, dpi="print")




# Function to read BED file and add sample name
read_bed <- function(file, sample_name) {
  data <- read_tsv(file, col_names = c("chromosome", "start", "end", "depth"))
  data$sample <- sample_name
  return(data)
}

# bars <- c("coat", "soap", "long", "sanger")
ana <- "/Users/kalebzuckerman/Desktop/paper_pushing/combined_analysis/longs"
file_suff <- "_1GB_DRCov.bed"

# Read in the BED files for four different samples
sample1 <- read_bed(paste(ana, "/b", "coat", file_suff, sep=""), "COATswga")
sample2 <- read_bed(paste(ana, "/b", "soap", file_suff, sep=""), "swga2.0")
sample3 <- read_bed(paste(ana, "/b", "sanger", file_suff, sep=""), "Sanger")

# Combine all datasets
# coverage_data <- bind_rows(sample1, sample2, sample3, sample4)
coverage_data <- bind_rows(sample1, sample2, sample3)

# Convert chromosome column to factor for proper ordering
coverage_data$chromosome <- factor(coverage_data$chromosome, levels = unique(coverage_data$chromosome))

# Read gene annotation file
genes <- read_tsv("/Users/kalebzuckerman/Desktop/paper_pushing/3D7_DR_genes.bed", col_names = c("chromosome", "gene_start", "gene_end", "gene"))

coverage_data <- coverage_data %>%
  left_join(genes %>% select(chromosome, gene))

# Plot using ggplot2
ridges <- ggplot(coverage_data, aes(xmin = start, xmax = end, ymin = 0, ymax = depth, fill = sample)) +
  geom_rect(size = 0.1) +
  # facet_grid(sample ~ chromosome, scales = "free_x", space = "free_x") +
  facet_grid(sample ~ gene, scales = "free", space = "free_x") +
  theme_minimal(base_size = 24) +
  labs(x = "Genomic Position",
       y = "Depth of Coverage",
       tag = "C") +
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
        plot.tag = element_text(size=36),
        # axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.margin = margin(t= 5.5, b=15, r=5.5, l=5.5),
        panel.spacing.x = unit(2, "lines"),
        panel.spacing.y = unit(1, "lines")
        # axis.title = element_text(size = 12)
  )

boxes <- ggplot(coverage_data, aes(x = sample, y = depth, fill = sample, weight=end-start)) +
  geom_boxplot(outlier.shape = NA) + # Boxplots with transparency
  scale_fill_manual(values = cs) +
  labs(y = "Depth of Coverage", fill = "Primer Set", tag="D") +
  facet_wrap(~gene, scales = "free", nrow = 1) +
  theme_minimal(base_size = 24) + 
  theme(axis.title.x = element_blank(),
        plot.tag = element_text(size=36),
        # axis.text = element_text(size = 12),
        strip.text = element_text(face = "bold"),
        plot.margin = margin(t=0, b=30, r=5.5, l=5.5),
        axis.text.x = element_blank(),
        legend.key.spacing.y = unit(0.5, "lines")
        # axis.title.y = element_text(size = 12)
  )

g <- grid.arrange(fig_2A + theme(plot.margin = margin(r = 40, l = 5.5, t=5.5, b=5.5)), fig_2B, ncol=2, widths=c(2,5))

sav <- grid.arrange(ridges, boxes, nrow=3, heights = c(5,2, 0.25))

final <- grid.arrange(g, sav, nrow=2)


ggsave("/Users/kalebzuckerman/Desktop/paper_pushing/final_figs/fig2_combined_final.jpg", plot = final, height = 6500, width=6000, units = "px")
# ggsave(paste0("/Users/kalebzuckerman/Desktop/placement_analysis/", "fig_DR.jpg"), plot = ridges,  width = 10, height = 7, )

