library(ggplot2)
library(gridExtra)

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
  theme_bw(base_size = 20) +
  theme(legend.position="none", plot.tag = element_text(size=22),
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
  geom_line() +
  theme_bw(base_size = 20) + 
  labs(tag="B", color="Primer Set") +
  scale_color_manual(values = cs) +
  scale_x_continuous(name ="Depth of Coverage", breaks=c(1, 5, seq(10,50, by=10)), expand=c(.01,0), limits=c(0,51)) +
  scale_y_continuous( name = "Percent of Genome Covered", breaks=seq(0,100, by = 20), limits = c(0,100)) +
  geom_vline(xintercept = 10, linetype = "dashed", color="grey") + 
  geom_vline(xintercept = 5, linetype = "dashed", color="grey") + 
  theme(plot.tag=element_text(size=22))
fig_2B

g <- grid.arrange(fig_2A + theme(plot.margin = margin(r = 40, l = 5.5, t=5.5, b=5.5)), fig_2B, ncol=2, widths=c(2,5))
ggsave(file=paste("/Users/kalebzuckerman/Desktop/paper_pushing/final_figs","/fig2.jpeg", sep=""), plot = g, units="px", width=4400, height=2120, dpi="print")



