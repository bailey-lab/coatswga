library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library(ggstats)
library(gridExtra)

# cs <- c("#C8102EFF", "#69B3E7FF", "#0E3492FF")#, "#FFC72CFF")
cs <- c("COATswga" = "#EE6100FF", "swga2.0" = "#1BB6AFFF", "Sanger" = "#9093A2FF")

df <- data.frame(matrix(ncol=5, nrow=0))

bars <- c("bcoat", "bsoap", "bsanger")
ana <- "/Users/kalebzuckerman/Desktop/paper_pushing/combined_analysis"
file_suff <- "_1GB_bin10.sorted.bedgraph"


ord <- paste0("chr", seq(1,14))
names <- c("COATswga", "swga2.0", "Sanger")
for (b in bars) {
  file_path <- paste(ana, "/", b, file_suff, sep="")
  f <- read.table(file_path)
  f <- f[!(f$V1 %in% c('chrM',"chrP")),]
  f$V5 <- names[which(bars == b)]
  df <- rbind(df, f)
}
df <- df %>% rename("chr" = V1, "start" = V2, "end" = V3, "depth" = V4, "sample" = V5)
df$chr <- factor(df$chr, levels = ord)
df$sample <- factor(df$sample, levels = names)
df <- df[order(df$chr, df$sample),]

# df <- df %>% filter(chr %in% c("chr1"))

fplot1 <- ggplot(df) +
  geom_rect(
    aes(xmin=start/1000000,xmax=end/1000000,ymin=0,ymax=depth, fill=sample), 
    linewidth = 0.1
    ) +
  coord_cartesian(ylim=c(0,50)) +
  theme_bw(base_size = 18) +
  scale_y_continuous(name="Depth of Coverage",expand=c(0,0)) +
  # scale_color_manual(values = cs, ) +
  scale_fill_manual(values = cs) +
  labs(tag = "A") +
  facet_grid(sample ~ chr, scales="free_x", space="free_x") +
  theme(panel.spacing.x = unit(0, "lines"), 
        panel.spacing.y = unit(0.75, "lines"), 
        legend.position = "none",
        axis.text = element_text(size = 12),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(face = "bold")) +
  guides(x="none")
# fplot1

b <- ggplot(df, aes(x = sample, y = depth, fill = sample, weight=(end - start)/10)) +
  geom_boxplot(outlier.shape = NA) + # Boxplots with transparency
  facet_wrap(~chr, scales = "free", nrow = 1) +
  coord_cartesian(ylim=c(0,100)) +
  theme_minimal(base_size = 18) +
  scale_y_continuous(name="Depth of Coverage",expand=c(0,0)) +
  labs(tag = "B", fill = "Primer Set") +
  # scale_color_manual(values = cs) +
  scale_fill_manual(values = cs) +
  theme(panel.spacing.x = unit(0, "lines"), 
        panel.spacing.y = unit(0.75, "lines"), 
        # legend.position = "none",
        # strip.background = element_rect(fill = "white"),
        strip.text = element_text(face = "bold"),
        axis.text = element_text(size = 12),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        plot.margin = margin(5.5, 5.5, 20, 5.5, "pt")) 
# b

arred <- grid.arrange(fplot1, b, nrow=2, heights=c(3,2))
ggsave(file="/Users/kalebzuckerman/Desktop/paper_pushing/final_figs/fig1_full.jpg", plot = arred, units="px", width=5000, height=2600, dpi="print")


