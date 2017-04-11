source("lib/functions.R")
source("lib/single-chr.R")

# load compiled C routines.
dyn.load("c-routines/do_exp_sums.so")
dyn.load("c-routines/calc_rel_dist.so")

# Read bed file of conserved regions (in this case "exons")
#bed = read.bed("primateElements.X")
gbed = read.bed.gmap("yri.piD.posInCm")

chrom = 'chr22';
args <- commandArgs(T)
par = c(0.768, 0.0001, 150000)

if (! is.na(args[3])) par <- as.numeric(args)
bed <- gbed[ gbed$chrom == chrom,]

# simulate some diversity data using our model (and add some noise)
# par is a vector of parameter values = A (neutral diversity), B (reduction per element) 
# and d (exponential parameter controlling scale over which diversity recovers).


# generate some positions along the chromosome
pos <- seq( head(bed$gStart, n=1) + 1000000, tail(bed$gEnd,n=1)-1000000, 1000)

# and make some fake data
sim = get_expected_gmap(pos, A=par[1],  B=par[2], d=par[3], bed=bed, max=1e6, mult=T, trunc=T)
# make some noise
sim <- sim / rnorm(length(sim), mean=1, sd=0.1)

# and use the longer-running time way to infer par
a <- quantile(sim, 0.99)
b <- quantile(sim, 0.01)


effElemLength<-1000

p3Min <- 2^5
p3Max <- 2^30

stop <- 10

results <- data.frame(matrix(ncol = 4, nrow = stop))
colnames(results) <- c('sumsq', 'p1', 'p2', 'p3')
x <- 1

p1 <- a
p2Start <- log(b/p1)/ -effElemLength


for (i in 0:stop) {
  p1 <- a + a*i/20
  p2 <- log(b/p1)/ -effElemLength
  
  
  best <- do_optimisation_gmap(pos, sim, bed, A=p1, B=p2) 
  print(best$par)
  
  best2<- optim(best$par, function(x) { sum_sq_gmap(pos, sim, A=x[1], B=x[2], d=x[3], bed=bed, max=1e6,mult=T) } ) 
 
  results[x,] <- c(best2$value, best2$par[1], best2$par[2], best2$par[3])
  print(results[x,])
  
  x<- x + 1   
}


