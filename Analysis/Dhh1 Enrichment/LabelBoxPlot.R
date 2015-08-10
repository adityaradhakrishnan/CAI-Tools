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

Dat      <- read.table('Test.txt',header=TRUE)

WTRO     <- Dat$WTProf/Dat$WTRNA
KORO     <- Dat$KOProf/Dat$KORNA
OERO     <- Dat$OEProf/Dat$OERNA
CAI      <- Dat$CAI

KOEnhance <- KORO/WTRO
OEEnhance <- OERO/WTRO

ROFrame <- data.frame(CAI,WTRO,KORO,OERO)
EnFrame <- data.frame(CAI,KOEnhance,OEEnhance)

Dat.RNA  <- melt(Dat, id.vars = 'CAI', measure.vars=c('WTRNA','KORNA','OERNA'))
Dat.Prof <- melt(Dat, id.vars = 'CAI', measure.vars=c('WTProf','KOProf','OEProf'))
Dat.RO   <- melt(ROFrame, id.vars = 'CAI', measure.vars=c('WTRO','KORO','OERO'))
Dat.En   <- melt(EnFrame, id.vars = 'CAI', measure.vars=c('KOEnhance','OEEnhance'))

ggplot(EnFrame, aes(x = factor(CAI), y = OEEnhance)) + geom_violin() + scale_y_log10(limits = c(0.075,15)) + theme_bw() + labs(x = 'CAI', y = 'Ribosome Occupancy (Dhh1 OE Relative to WT)') + stat_summary(fun.data = give.n, geom = "text", fun.y = median, vjust = 0.5, size=4) + theme(text = element_text(size=16))

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