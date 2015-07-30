# function for number of observations 
give.n <- function(x){
  return(c(y = median(x)*1.05, label = length(x))) 
  # experiment with the multiplier to find the perfect position
}

# function for mean labels
mean.n <- function(x){
  return(c(y = median(x)*0.97, label = round(mean(x),2))) 
  # experiment with the multiplier to find the perfect position
}

ggplot(a, aes(factor(V2), V4/V5)) + geom_boxplot() + scale_y_log10(limits = c(0.01,100)) + stat_summary(fun.data = give.n, geom = "text", fun.y = median, vjust = -0.5, size = 3) + stat_summary(fun.data = mean.n, geom='text', fun.y = mean, colour = 'red', vjust = 1.5, size = 3) + labs(x='CAI', y = 'Enrichment in Ribosome Footprings in OE/WT')