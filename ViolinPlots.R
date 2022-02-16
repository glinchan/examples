#!/usr/bin/env Rscript


##
# IMPORT LIBRARIES
##

suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(scales))
suppressPackageStartupMessages(require(gridExtra))
suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(dplyr))
suppressPackageStartupMessages(require(ggpubr))

##
# DEFINE ARGUMENTS
##

option_list = list(
  make_option(c("-c", "--csv"), action="store", type="character", default="data.csv", help="Table (in comma-separated values) [default %default]", metavar="string"),
  make_option(c("-d", "--csv2"), action="store", type="character", default ="data.csv",help="Table (in comma-separated values) [default %default]", metavar ="string"));

opt = parse_args(OptionParser(option_list=option_list))
fileName  <- opt$csv
fileName2 <- opt$csv2


##
# READ DATA
##

B <- read.csv(fileName)
C <- read.csv(fileName2)

## Violin Plots

compareCON<-unique(B$Consensus_ID)

s <- ggplot(B, aes(x=Consensus_ID, y=Substitutions, fill=Comparison)) + geom_violin()+theme(legend.position="none")+theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.title.y = element_blank())

s <- s + ggtitle("A")
s <- s + geom_boxplot(width=0.1, fill = "white",outlier.shape=NA)
s <- s + stat_compare_means(tip.length = 0,method = "wilcox.test",comparisons = list(c(compareCON[1],compareCON[2]),c(compareCON[3],compareCON[4]),c(compareCON[5],compareCON[6]),c(compareCON[7],compareCON[8]),c(compareCON[9],compareCON[10]),c(compareCON[11],compareCON[12]),c(compareCON[13],compareCON[14]),c(compareCON[15],compareCON[16])), paired = FALSE,label.y=500,label = "p.format")
s <- s + scale_fill_manual(values = c("#FF6666", "#3399FF"))
s <- s + theme(text = element_text(size=14))
compareCON2<-unique(C$Consensus_ID)
s2 <- ggplot(C, aes(x=Consensus_ID, y=Substitutions, fill=Comparison)) + geom_violin()+theme(legend.position="none")+theme(axis.title.x=element_blank(),axis.title.y = element_blank())

s2 <- s2 + ggtitle("B")
s2 <- s2 + geom_boxplot(width=0.1, fill = "white",outlier.shape=NA)
s2 <- s2 + stat_compare_means(tip.length = 0, method = "wilcox.test",comparisons = list(c(compareCON2[1],compareCON2[2]),c(compareCON2[3],compareCON2[4]),c(compareCON2[5],compareCON2[6]),c(compareCON2[7],compareCON2[8]),c(compareCON2[9],compareCON2[10]),c(compareCON2[11],compareCON2[12]),c(compareCON2[13],compareCON2[14]),c(compareCON2[15],compareCON2[16])), paired = FALSE,label.y=500,label = "p.format")
s2 <- s2 + scale_fill_manual(values = c("#FFCC00", "#3399FF"))
s2 <- s2 + theme(text = element_text(size=14))
g <- grid.arrange(arrangeGrob(s, s2, nrow=2))
ggsave(filename=paste("Substitutions_Vplots_",fileName,".png",sep=""), g, width=20, height=15, dpi=400)

#deletions
d <- ggplot(B, aes(x=Consensus_ID, y=Deletions, fill=Comparison)) + geom_violin()+theme(legend.position="none")+theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.title.y = element_blank())
d <- d + ggtitle("A")
d <- d + geom_boxplot(width=0.1, fill = "white",outlier.shape=NA)
d <- d + stat_compare_means(tip.length = 0,method = "wilcox.test",comparisons = list(c(compareCON[1],compareCON[2]),c(compareCON[3],compareCON[4]),c(compareCON[5],compareCON[6]),c(compareCON[7],compareCON[8]),c(compareCON[9],compareCON[10]),c(compareCON[11],compareCON[12]),c(compareCON[13],compareCON[14]),c(compareCON[15],compareCON[16])), paired = FALSE,label.y=100,label = "p.format")
d <- d + scale_fill_manual(values = c("#FF6666", "#3399FF"))
d <- d + theme(text = element_text(size=14))
d2 <- ggplot(C, aes(x=Consensus_ID, y=Deletions, fill=Comparison)) + geom_violin()+theme(legend.position="none")+theme(axis.title.x=element_blank(),axis.title.y = element_blank())
d2 <- d2 + ggtitle("B")
d2 <- d2 + geom_boxplot(width=0.1, fill = "white",outlier.shape=NA)
d2 <- d2 + stat_compare_means(tip.length = 0, method = "wilcox.test",comparisons = list(c(compareCON2[1],compareCON2[2]),c(compareCON2[3],compareCON2[4]),c(compareCON2[5],compareCON2[6]),c(compareCON2[7],compareCON2[8]),c(compareCON2[9],compareCON2[10]),c(compareCON2[11],compareCON2[12]),c(compareCON2[13],compareCON2[14]),c(compareCON2[15],compareCON2[16])), paired = FALSE,label.y=100,label = "p.format")
d2 <- d2 + scale_fill_manual(values = c("#FFCC00", "#3399FF"))
d2 <- d2 + theme(text = element_text(size=14))
g2 <- grid.arrange(arrangeGrob(d,d2, nrow=2))
ggsave(filename=paste("Deletions_Vplots_",fileName,".png",sep=""), g2, width=20, height=15, dpi=400)

##insertions
i <- ggplot(B, aes(x=Consensus_ID, y=Insertions, fill=Comparison)) + geom_violin()+theme(legend.position="none")+theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.title.y = element_blank())
i <- i + ggtitle("A")
i <- i + geom_boxplot(width=0.1, fill = "white",outlier.shape=NA)
i <- i + stat_compare_means(tip.length = 0,method = "wilcox.test",comparisons = list(c(compareCON[1],compareCON[2]),c(compareCON[3],compareCON[4]),c(compareCON[5],compareCON[6]),c(compareCON[7],compareCON[8]),c(compareCON[9],compareCON[10]),c(compareCON[11],compareCON[12]),c(compareCON[13],compareCON[14]),c(compareCON[15],compareCON[16])), paired = FALSE,label.y=100,label = "p.format")
i <- i + scale_fill_manual(values = c("#FF6666", "#3399FF"))
i <- i + theme(text = element_text(size=14))

i2 <- ggplot(C, aes(x=Consensus_ID, y=Insertions, fill=Comparison)) + geom_violin()+theme(legend.position="none")+theme(axis.title.x=element_blank(),axis.title.y = element_blank())
i2 <- i2 + ggtitle("B")
i2 <- i2 + geom_boxplot(width=0.1, fill = "white",outlier.shape=NA)
i2 <- i2 + stat_compare_means(tip.length = 0, method = "wilcox.test",comparisons = list(c(compareCON2[1],compareCON2[2]),c(compareCON2[3],compareCON2[4]),c(compareCON2[5],compareCON2[6]),c(compareCON2[7],compareCON2[8]),c(compareCON2[9],compareCON2[10]),c(compareCON2[11],compareCON2[12]),c(compareCON2[13],compareCON2[14]),c(compareCON2[15],compareCON2[16])), paired = FALSE,label.y=100,label = "p.format")
i2 <- i2 + scale_fill_manual(values = c("#FFCC00", "#3399FF"))
i2 <- i2 + theme(text = element_text(size=14))
g3 <- grid.arrange(arrangeGrob(i,i2, nrow=2))
ggsave(filename=paste("Insertions_Vplots_",fileName,".png",sep=""), g3, width=20, height=15, dpi=400)



