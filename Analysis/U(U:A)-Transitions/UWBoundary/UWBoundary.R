library(ggplot2)
require(gridExtra)

lm_eqn = function(m) {

  l <- list(a = format(coef(m)[1], digits = 2),
      b = format(abs(coef(m)[2]), digits = 2),
      r2 = format(summary(m)$r.squared, digits = 3));

  if (coef(m)[2] >= 0)  {
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,l)
  } else {
    eq <- substitute(italic(y) == a - b %.% italic(x)*","~~italic(r)^2~"="~r2,l)    
  }

  as.character(as.expression(eq));                 
}

Protein <- read.table('UWBoundary-Protein.tsv', header=TRUE)

CountFit <- lm(Expression ~ UWCounts, Protein)
FracFit <- lm(Expression ~ UWFrac, Protein)

plotCounts <- ggplot() + geom_point(data = Protein, aes(x = UWCounts, y = Expression), cex = 1) + theme_bw() + labs(x = 'U(A/U) Codon Boundaries in Gene', y = 'Protein Expression (Molecules/Cell)') + scale_y_log10(limits = c(10,10000000)) + scale_x_log10(limits = c(1,1000)) + annotate("text", x = 7, y = 7000000, label = lm_eqn(CountFit), parse = TRUE, cex = 4)

plotFrac <- ggplot() + geom_point(data = Protein, aes(x = UWFrac, y = Expression), cex = 1) + theme_bw() + labs(x = 'Fraction of Codon Boundaries in Gene That Are U(A/U)', y = 'Protein Expression (Molecules/Cell)') + scale_y_log10(limits = c(10,10000000)) + annotate("text", x = .11, y = 7000000, label = lm_eqn(FracFit), parse = TRUE, cex = 4)

pdf('Protein.pdf', width=11, height=5)

grid.arrange(plotCounts, plotFrac, ncol=2)

dev.off()

RNA <- read.table('UWBoundary-RNA.tsv', header=TRUE)

PolyACountFit <- lm(PolyAHL ~ UWCounts, RNA)
PolyAFracFit  <- lm(PolyAHL ~ UWFrac, RNA)
CountFit      <- lm(HL ~ UWCounts, RNA)
FracFit       <- lm(HL ~ UWFrac, RNA)

plotPolyACounts <- ggplot() + geom_point(data = RNA, aes(x = UWCounts, y = PolyAHL), cex = 1) + theme_bw() + labs(x = 'U(A/U) Codon Boundaries in Gene', y = 'Poly(A) Half Life (Minutes)') + scale_y_log10(limits = c(1,60)) + scale_x_log10(limits = c(1,1000)) + annotate("text", x = 7, y = 55, label = lm_eqn(PolyACountFit), parse = TRUE, cex = 4)

plotPolyAFrac <- ggplot() + geom_point(data = RNA, aes(x = UWFrac, y = PolyAHL), cex = 1) + theme_bw() + labs(x = 'Fraction of Codon Boundaries in Gene That Are U(A/U)', y = 'Poly(A) Half Life (Minutes)') + scale_y_log10(limits = c(1,60))  + annotate("text", x = .11, y = 55, label = lm_eqn(PolyAFracFit), parse = TRUE, cex = 4)

plotCounts <- ggplot() + geom_point(data = RNA, aes(x = UWCounts, y = HL), cex = 1) + theme_bw() + labs(x = 'U(A/U) Codon Boundaries in Gene', y = 'Half Life (Minutes)')  + scale_y_log10(limits = c(1,80)) + scale_x_log10(limits = c(1,1000)) + annotate("text", x = 5.4, y = 75, label = lm_eqn(CountFit), parse = TRUE, cex = 4)

plotFrac <- ggplot() + geom_point(data = RNA, aes(x = UWFrac, y = HL), cex = 1) + theme_bw() + labs(x = 'Fraction of Codon Boundaries in Gene That Are U(A/U)', y = 'Half Life (Minutes)') + scale_y_log10(limits = c(1,80))  + annotate("text", x = .11, y = 75, label = lm_eqn(FracFit), parse = TRUE, cex = 4)

pdf('RNA.pdf', width=11, height=11)

grid.arrange(plotPolyACounts, plotPolyAFrac, plotCounts, plotFrac, ncol=2, nrow=2)

dev.off()
