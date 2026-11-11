library(ggplot2)
library(dplyr)
library(gridExtra)
library(fuzzyjoin)
library(plyr)
library(readr)

cs <- c("COATswga" = "#EE6100FF", "swga2.0" = "#1BB6AFFF", "Sanger" = "#9093A2FF")


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

boxes <- ggplot(coverage_data, aes(x = sample, y = depth, fill = sample, weight=end-start)) +
  geom_boxplot(outlier.shape = NA) + # Boxplots with transparency
  scale_fill_manual(values = cs) +
  labs(y = "Depth of Coverage", fill = "Primer Set", tag="B") +
  facet_wrap(~gene, scales = "free", nrow = 1) +
  theme_minimal(base_size = 20) + 
  theme(axis.title.x = element_blank(),
        axis.text = element_text(size = 12),
        strip.text = element_text(face = "bold"),
        plot.margin = margin(t=0, b=30, r=5.5, l=5.5),
        axis.text.x = element_blank(),
        # axis.title.y = element_text(size = 12)
        )
library(dplyr)

dr <- read_csv("/Users/kalebzuckerman/Desktop/paper_pushing/keydrugresistancegenesandloci.csv") %>% 
  select(CHROM, POS, mutation_name)  %>%
  unique()

combined <- inner_join(dr, coverage_data, by = c("CHROM" = "chromosome"), relationship = "many-to-many") %>%
  filter(POS >= start, POS < end) %>%
  select(CHROM, POS, mutation_name, depth, sample)

maxed <- combined %>% 
  dplyr::group_by(mutation_name) %>%
  dplyr::summarise(sample_max = sample[[which.max(depth)]])

max_dat <- data.frame(table(maxed$sample_max)) 
# %>%
# add_row(Var1 = c("COATxSanger", "COATxSOAP"), Freq = c(0,0))

br <- ggplot(max_dat, aes(x = Var1, y = Freq, fill = Var1)) + 
  geom_bar(stat = "identity") +
  labs(y = expression("   # of Loci With\nHighest Coverage"), tag = "C") +
  scale_fill_manual(values = cs) + 
  theme_minimal(base_size = 20) + 
  theme(legend.position = "none", 
        axis.text = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16)
        )
br
sav <- grid.arrange(ridges, boxes,  arrangeGrob(br, ncol=2), nrow=3, heights = c(5,2, 2))

ggsave("/Users/kalebzuckerman/Desktop/paper_pushing/final_figs/figDR.jpg", plot = sav, height = 3500, width=4000, units = "px")
# ggsave(paste0("/Users/kalebzuckerman/Desktop/placement_analysis/", "fig_DR.jpg"), plot = ridges,  width = 10, height = 7, )

