library(ggplot2)
library(gridExtra)

cs <- c("COATswga" = "#EE6100FF", "swga2.0" = "#1BB6AFFF", "Sanger" = "#9093A2FF")


bars <- c("coat", "soap", "sanger")
ana <- "/Users/kalebzuckerman/Desktop/paper_pushing/combined_analysis/stepped"
names <- c("COATswga", "swga2.0", "Sanger")


dat_1X <- data.frame(matrix(ncol=4, nrow=0))

for (b in bars){
  file_path <- paste(ana, "/", b, "_1X_stepped_cov.csv", sep="")
  f <- read.csv(file_path, header=FALSE)
  f$V3 <- names[which(bars == b)]
  dat_1X <- rbind(dat_1X, f)
}
dat_1X

steps_1X <- ggplot(dat_1X) +
  geom_line(aes(x=V1, y=V2, color=V3), size=0.75) +
  # geom_smooth(aes(x=V1, y=V2, color=V3), method = "loess", se = FALSE, span=0.5, size=0.5) +
  scale_x_continuous(name ="Sequencing effort (Mbp sequenced)", expand=c(.01,0), limits=c(10,1010), breaks=c(10, 50, seq(100,1000, by=100) )) +
  scale_y_continuous( name = "Percent of Genome Covered (1X)", breaks=seq(0,100, by = 20), limits = c(0,100)) +
  scale_color_manual(values = cs) +
  labs(color="Primer Set", tag="A") +
  theme_bw(base_size = 20) + 
  theme(panel.grid.minor= element_blank(), 
        legend.position = "none", 
        # axis.title = element_text(size = 14),
        axis.text = element_text(size = 12, angle=60, hjust=1),
        plot.tag = element_text(size=22))

steps_1X  


dat_5X <- data.frame(matrix(ncol=4, nrow=0))

for (b in bars){
  file_path <- paste(ana, "/", b, "_5X_stepped_cov.csv", sep="")
  f <- read.csv(file_path, header=FALSE)
  f$V3 <- names[which(bars == b)]
  dat_5X <- rbind(dat_5X, f)
}
dat_5X


steps_5X <- ggplot(dat_5X) +
  geom_line(aes(x=V1, y=V2, color=V3), size=0.75) +
  # geom_smooth(aes(x=V1, y=V2, color=V3), method = "loess", se = FALSE, span=0.5, size=0.5) +
  # facet_wrap(~V3) +
  # scale_x_continuous(name ="Sequencing effort (Mbp sequenced)", expand=c(.01,0), limits=c(10,500), breaks=c(10, 50, seq(100,600, by=100) )) +
  scale_x_continuous(name ="Sequencing effort (Mbp sequenced)", expand=c(.01,0), limits=c(10,1010), breaks=c(10, 50, seq(100,1000, by=100) )) +
  scale_y_continuous( name = "Percent of Genome Covered (5X)", breaks=seq(0,100, by = 20), limits = c(0,100)) +
  scale_color_manual(values = cs) +
  labs(color="Primer Set", tag="B") +
  theme_bw(base_size = 20) + 
  theme(panel.grid.minor= element_blank(), 
        legend.position = "none", 
        # axis.title = element_text(size = 14),
        axis.text = element_text(size = 12, angle=60, hjust=1),
        plot.tag = element_text(size=22))

steps_5X  

dat_10X <- data.frame(matrix(ncol=3, nrow=0))
for (b in bars){
  file_path <- paste(ana, "/", b, "_10X_stepped_cov.csv", sep="")
  f <- read.csv(file_path, header=FALSE)
  f$V3 <- names[which(bars == b)]
  dat_10X <- rbind(dat_10X, f)
}
dat_10X

steps_10X <- ggplot(dat_10X, aes(x=V1, y=V2, color=V3)) +
  geom_line(size=0.75) +
  # scale_x_continuous(name ="Sequencing effort (Mbp sequenced)", expand=c(.01,0), limits=c(10,610), breaks=c(10, 50, seq(100,600, by=100) )) +
  scale_x_continuous(name ="Sequencing effort (Mbp sequenced)", expand=c(.01,0), limits=c(10,1010), breaks=c(10, 50, seq(100,1000, by=100) )) +
  scale_y_continuous( name = "Percent of Genome Covered (10X)", breaks=seq(0,100, by = 20), limits = c(0,100)) +
  scale_color_manual(values = cs) +
  labs(color="Primer Set", tag="C") +
  theme_bw(base_size = 20) + 
  theme(panel.grid.minor= element_blank(), 
        # axis.title = element_text(size = 14),
        axis.text = element_text(size = 12, angle=60, hjust=1),
        plot.tag = element_text(size=22))

steps_10X  

fig3 <- grid.arrange(steps_1X, steps_5X, steps_10X, ncol=3, widths=c(4, 4, 5))
fig3
ggsave(file="/Users/kalebzuckerman/Desktop/paper_pushing/final_figs/fig_stepped.jpg", plot = fig3, units="px", width=6400, height=2000, dpi="print")


