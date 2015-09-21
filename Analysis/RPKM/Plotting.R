Un <- read.table('Unmod.txt')
IF <- read.table('Iter-FirstMod.txt')
IS <- read.table('Iter-SecondMod.txt')

pdf('CDF-Unmodified.pdf')

ggplot(data = Un) + stat_ecdf(aes(x = 1e9*V6/(V2*sum(V6)))) + stat_ecdf(aes(x = 1e9*V7/(V2*sum(V7))), colour = 'blue') + stat_ecdf(aes(x = 1e9*V8/(V2*sum(V8))), colour = 'red') + scale_x_log10() + theme_bw() + xlab('mRNA Seq (RPKM)') + ylab('Cumulative Fraction')

dev.off()

pdf('CDF-Iter-First.pdf')

ggplot(data = IF) + stat_ecdf(aes(x = 1e9*V6/(V2*sum(V6)))) + stat_ecdf(aes(x = 1e9*V7/(V2*sum(V7))), colour = 'blue') + stat_ecdf(aes(x = 1e9*V8/(V2*sum(V8))), colour = 'red') + scale_x_log10() + theme_bw() + xlab('mRNA Seq (RPKM)') + ylab('Cumulative Fraction')

dev.off()

pdf('CDF-Iter-Second.pdf')

ggplot(data = IS) + stat_ecdf(aes(x = 1e9*V6/(V2*sum(V6)))) + stat_ecdf(aes(x = 1e9*V7/(V2*sum(V7))), colour = 'blue') + stat_ecdf(aes(x = 1e9*V8/(V2*sum(V8))), colour = 'red') + scale_x_log10() + theme_bw() + xlab('mRNA Seq (RPKM)') + ylab('Cumulative Fraction')

dev.off()

WTRPKM <- 1e9*Un$V6/(Un$V2*sum(Un$V6))
KORPKM <- 1e9*Un$V7/(Un$V2*sum(Un$V7))
OERPKM <- 1e9*Un$V8/(Un$V2*sum(Un$V8))

DensityCheck <- (WTRPKM > 10)*(KORPKM > 10)*(OERPKM > 10)
mode(DensityCheck) <- "logical"

KO  <- KORPKM[DensityCheck]
WT  <- WTRPKM[DensityCheck]
OE  <- OERPKM[DensityCheck]
CAI <- Un$V3[DensityCheck]

length(CAI[CAI == 0.3])
length(CAI[CAI == 0.4])
length(CAI[CAI == 0.5])
length(CAI[CAI == 0.6])
length(CAI[CAI == 0.7])
length(CAI[CAI == 0.8])
length(CAI[CAI == 0.9])

ggplot() + geom_violin(aes(factor(CAI),KO/WT)) + scale_y_log10(limits = c(0.1,10)) + theme_bw()

ggplot() + geom_point(aes(x = KO/WT, y = OE/WT)) + scale_x_log10() + scale_y_log10() + theme_bw()

WTRPKM <- 1e9*IF$V6/(IF$V2*sum(IF$V6))
KORPKM <- 1e9*IF$V7/(IF$V2*sum(IF$V7))
OERPKM <- 1e9*IF$V8/(IF$V2*sum(IF$V8))

KO <- KORPKM[WTRPKM > 10]
WT <- WTRPKM[WTRPKM > 10]
OE <- OERPKM[WTRPKM > 10]
CAI <- IF$V3[WTRPKM > 10]

ggplot() + geom_point(aes(x = KO/WT, y = OE/WT)) + scale_x_log10() + scale_y_log10() + theme_bw()
ggplot() + geom_violin(aes(factor(CAI),OE/WT)) + scale_y_log10(limits = c(0.1,20)) + theme_bw()

WTRPKM <- 1e9*IS$V6/(IS$V2*sum(IS$V6))
KORPKM <- 1e9*IS$V7/(IS$V2*sum(IS$V7))
OERPKM <- 1e9*IS$V8/(IS$V2*sum(IS$V8))

KO <- KORPKM[WTRPKM > 10]
WT <- WTRPKM[WTRPKM > 10]
OE <- OERPKM[WTRPKM > 10]
CAI <- IS$V3[WTRPKM > 10]

ggplot() + geom_point(aes(x = KO/WT, y = OE/WT)) + scale_x_log10() + scale_y_log10()
ggplot() + geom_violin(aes(factor(CAI),OE/WT)) + scale_y_log10(limits = c(0.1,20))
