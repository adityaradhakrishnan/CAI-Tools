# function for number of observations 
give.n <- function(x){
  return(c(y = median(x)*0.97, label = length(x))) 
  # experiment with the multiplier to find the perfect position
}

# function for mean labels
mean.n <- function(x){
  return(c(y = median(x)*0.97, label = round(mean(x),2))) 
  # experiment with the multiplier to find the perfect position
}

library(ggplot2)

Dat      <- read.table('RPKM.txt')

CAI <- Dat$V2
Dh1 <- Dat$V3
Dh2 <- Dat$V4
Sb1 <- Dat$V5
Sb2 <- Dat$V6
Ls1 <- Dat$V7
Ls2 <- Dat$V8
Pt1 <- Dat$V9
Pt2 <- Dat$V10
WT  <- Dat$V11

Frame <- data.frame(CAI,Dh1,Dh2,Sb1,Sb2,Ls1,Ls2,Pt1,Pt2,WT)

ggplot() + geom_violin(aes(x = factor(CAI), y = Dh1/WT), colour = 'blue', alpha = 0.5) + geom_violin(aes(x = factor(CAI), y = Dh2/WT), colour = 'orange', alpha = 0.5) + scale_y_log10(limits = c(0.01,100)) + theme_bw() + labs(x = 'CAI', y = 'mRNA Enrichment in Dhh1 Clip-SEQ Relative to Wild Type')

ggplot(aes(x = factor(CAI), y = Dh1)) + geom_violin() + scale_y_log10(limits = c(0.075,15)) + theme_bw() + labs(x = 'CAI', y = 'Ribosome Occupancy (Dhh1 OE Relative to WT)') + stat_summary(fun.data = give.n, geom = "text", fun.y = median, vjust = 0.5, size=4) + theme(text = element_text(size=16))

ggplot(Dat, aes(x = factor(CAI), y = Dat$KORNA/Dat$WTRNA)) + geom_violin() + scale_y_log10(limits = c(0.075,15)) + theme_bw() + labs(x = 'CAI', y = 'Steady State mRNA Levels (Dhh1 KO Relative to WT)') + stat_summary(fun.data = give.n, geom = "text", fun.y = median, vjust = 0.5, size=4) + theme(text = element_text(size=16))

ggplot(Dat, aes(x = factor(CAI), y = Dat$KORNA/Dat$WTRNA*Dat$CAI^0.35)) + geom_violin() + scale_y_log10(limits = c(0.075,15)) + theme_bw() + labs(x = 'CAI', y = 'Steady State mRNA Levels (Dhh1 OE Relative to WT)') + stat_summary(fun.data = give.n, geom = "text", fun.y = median, vjust = 0.5, size=4) + theme(text = element_text(size=16))

ggplot() + geom_point(dat = EnFrame, aes(x = KOREN, y = OEEnhance), cex = 1) + scale_y_log10(limits = c(0.075,15)) + scale_x_log10(limits = c(0.05,30)) + theme_bw() + labs(x = 'Steady State mRNA Levels (Dhh1 KO Relative to WT)', y = 'Ribosome Occupancy (Dhh1 OE Relative to WT)') + theme(text = element_text(size=16)) + geom_point(dat = ModFrame, aes(x = KOREN3, y = OEEnhance3), cex = 2, colour = 'blue')  + geom_point(dat = ModFrame2, aes(x = KOREN9, y = OEEnhance9), cex = 2, colour = 'orange') 

ggplot(EnFrame, aes(x = factor(CAI), y = KOEnhance)) + geom_violin() + scale_y_log10(limits = c(0.075,15)) + theme_bw() + labs(x = 'CAI', y = 'Ribosome Occupancy (Dhh1 KO Relative to WT)') + stat_summary(fun.data = give.n, geom = "text", fun.y = median, vjust = 0.5, size=4) + theme(text = element_text(size=16))

ggplot(EnFrame, aes(x = factor(CAI), y = OEEnhance)) + geom_boxplot() + scale_y_log10(limits = c(0.075,15)) + theme_bw() + labs(x = 'CAI', y = 'Ribosome Occupancy (Dhh1 OE Relative to WT)') + stat_summary(fun.data = give.n, geom = "text", fun.y = median, vjust = -0.4, size=3) + theme(text = element_text(size=16))

ggplot(EnFrame, aes(x = factor(CAI), y = KOEnhance)) + geom_boxplot() + scale_y_log10(limits = c(0.075,15)) + theme_bw() + labs(x = 'CAI', y = 'Ribosome Occupancy (Dhh1 KO Relative to WT)') + stat_summary(fun.data = give.n, geom = "text", fun.y = median, vjust = 1.4, size=3) + theme(text = element_text(size=16))

ggplot(EnFrame, aes(x = factor(CAI), y = OEEnhance)) + geom_boxplot(outlier.size = 1) + scale_y_log10(limits = c(0.075,15)) + theme_bw() + labs(x = 'CAI', y = 'Ribosome Occupancy (Dhh1 OE Relative to WT)') + stat_summary(fun.data = give.n, geom = "text", fun.y = median, vjust = -0.4, size=3) + theme(text = element_text(size=16))

ggplot(EnFrame, aes(x = factor(CAI), y = KOEnhance)) + geom_boxplot(outlier.size = 1) + scale_y_log10(limits = c(0.075,15)) + theme_bw() + labs(x = 'CAI', y = 'Ribosome Occupancy (Dhh1 KO Relative to WT)') + stat_summary(fun.data = give.n, geom = "text", fun.y = median, vjust = 1.4, size=3) + theme(text = element_text(size=16))

ggplot(EnFrame, aes(x = factor(CAI), y = KOEnhance)) + geom_boxplot() + scale_y_log10(limits = c(0.075,15)) + theme_bw() + labs(x = 'CAI', y = 'Ribosome Occupancy (Dhh1 KO Relative to WT)') + stat_summary(fun.data = give.n, geom = "text", fun.y = median, vjust = 1.4, size=3) + theme(text = element_text(size=16))



ggplot(Dat.RNA) + geom_boxplot(aes(x = factor(CAI), y = value, color = variable)) + scale_y_log10() + theme_bw() + labs(x = 'CAI', y = 'mRNA Seq (RPKM)')
ggplot(Dat.Prof) + geom_boxplot(aes(x = factor(CAI), y = value, color = variable)) + scale_y_log10() + theme_bw() + labs(x = 'CAI', y = 'mRNA Seq (RPKM)')
ggplot(Dat.RO) + geom_boxplot(aes(x = factor(CAI), y = value, color = variable)) + scale_y_log10() + theme_bw() + labs(x = 'CAI', y = 'Ribosome Occupancy')
ggplot(EnFrame, aes(x = factor(CAI), y = KOEnhance)) + geom_boxplot(outlier.colour = 'grey80', outlier.size = 1) + scale_y_log10(limits = c(0.05,20)) + theme_bw() + labs(x = 'CAI', y = 'Ribosome Occupancy') + stat_summary(fun.data = give.n, geom = "text", fun.y = median, vjust = -1.5, size = 3) + stat_summary(fun.data = mean.n, geom='text', fun.y = mean, colour = 'red', vjust = 1.5, size = 3)



ggplot(a, aes(factor(V2), V4/V5)) + geom_boxplot() + scale_y_log10(limits = c(0.01,100)) + stat_summary(fun.data = give.n, geom = "text", fun.y = median, vjust = -0.5, size = 3) + stat_summary(fun.data = mean.n, geom='text', fun.y = mean, colour = 'red', vjust = 1.5, size = 3) + labs(x='CAI', y = 'Enrichment in Ribosome Footprings in OE/WT')