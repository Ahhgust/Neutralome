source("lib/functions.R")
source("lib/single-chr.R")

# load compiled C routines.
dyn.load("c-routines/do_exp_sums.so")
dyn.load("c-routines/calc_rel_dist.so")

# Read bed file of conserved regions (in this case "exons")
#bed = read.bed("primateElements.X")
bed = read.bed("chrom22.elements")
# simulate some diversity data using our model (and add some noise)
# par is a vector of parameter values = A (neutral diversity), B (reduction per element) 
# and d (exponential parameter controlling scale over which diversity recovers).
par = c(0.768, 0.0001, 150000)

par <- as.numeric(commandArgs(T))
if (is.na(par[3])) par = c(0.768, 0.0001, 150000) # defaults 

par <- par[1:3]

# generate some positions along the chromosome
pos <- seq( head(bed$chromStart, n=1) + 1000000, tail(bed$chromEnd,n=1)-1000000, 1000)

# and make some fake data
sim = get_expected(pos, A=par[1],  B=par[2], d=par[3], bed=bed, max=1e6, mult=T, trunc=T)
#sim = sim+rnorm(length(sim), mean=0, sd=0.1)
sim <- sim / rnorm(length(sim), mean=1, sd=0.1)
# and use the longer-running time way to infer par
a <- quantile(sim, 0.99)
b <- quantile(sim, 0.01)

#numSelectedSites <- sum( bed$chromEnd - bed$chromStart)

#expected <- 1/numSelectedSites*(log(sim)-log(a))

#distance <- pos
#for (i in 1:length(pos) )
#  distance[i] <- mean(abs(pos[i] - (bed$chromEnd + bed$chromStart)/2.0))

#fit <- lm(expected ~ distance)

#p1 <- a
#p2 <- coef(fit)[1]
#p3 <- -p2/coef(fit)[2]

#res = optim(c(a,b, par[3]), function(x) { sum_sq(pos, sim, A=x[1], B=x[2], d=x[3], bed=bed, max=1e6) } )
#res = do_optimisation(pos, sim, bed, A=a, B=b, mult=T,trunc=T)
#print( cat(par, par2sim, a, b))
#res

effElemLength<-1000

p3Min <- 2^5
p3Max <- 2^30

irange <- 21
jrange <- 20


df <- data.frame(matrix(ncol = 4, nrow = irange*jrange))
colnames(df) <- c('sumsq', 'p1', 'p2', 'p3')
x <- 1
fp <- function(x) { res = get_expected_gmap(pos, A=p1,  B=p2, d=x, bed=bed, max=1e6, mult=T, trunc=T); sum( (res-sim)^2) }

for (i in 0:20) {
  p1 <- a + a*i/20
  p2Start <- log(b/p1)/ -effElemLength

  for (j in seq(2,40,2)) {
    p2 <- p2Start / j
    best <- optimize(fp, c(p3Min, p3Max),tol=1)
    df[x,] <- c(best[[2]], p1, p2, best[[1]])
    print(df[x,])
    x <- x + 1
  }
}



