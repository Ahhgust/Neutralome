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
#par = c(0.768, 0.0001, 150000)

par <- as.numeric(commandArgs(T))
if (is.na(par[6])) stop("Oops! not enough arguments!")

par2sim <- par[4:6]
par <- par[1:3]

# generate some positions along the chromosome
pos <- seq( head(bed$chromStart, n=1) + 1000000, tail(bed$chromEnd,n=1)-1000000, 1000)

# and make some fake data
sim = get_expected(pos, A=par[1],  B=par[2], d=par[3], bed=bed, max=1e6, mult=T, trunc=T)

# and use the longer-running time way to infer par
a <- quantile(sim, 0.99)
b <- quantile(sim, 0.01)
res = optim(c(a,b, par[3]), function(x) { sum_sq(pos, sim, A=x[1], B=x[2], d=x[3], bed=bed, max=1e6) } )
#res = do_optimisation(pos, sim, bed, A=a, B=b, mult=T,trunc=T)
print( cat(par, par2sim, a, b))
res





