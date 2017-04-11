source("lib/functions.R")
source("lib/single-chr.R")

# load compiled C routines.
dyn.load("c-routines/do_exp_sums.so")
dyn.load("c-routines/calc_rel_dist.so")

# function to plot an conserved element on a figure
addexon = function(r, ylim, ...) {
   polygon(c(r[2],r[3],r[3],r[2]), c(ylim[1], ylim[1], ylim[2], ylim[2]), ...)
}

# Read bed file of conserved regions (in this case "exons")
bed = read.bed("sample/exons.bed")

# simulate some diversity data using our model (and add some noise)
# par is a vector of parameter values = A (neutral diversity), B (reduction per element) 
# and d (exponential parameter controlling scale over which diversity recovers).
par = c(0.768, 0.0001, 150000)

# pos is a vector of positions at which diversity is to be predicted.
pos = seq(2254500, 5089500, 1000)

# generate estimates of predicted diversity
sim = get_expected(pos, A=par[1],  B=par[2], d=par[3], bed=bed, max=1e6)

# add noise to diversity estimates
sim = sim+rnorm(length(sim), mean=0, sd=0.1)


# Now to estimate the parameters from the data by minimising the sum of squares
# do full optimisation over all three parameters
res = optim(par, function(x) { sum_sq(pos, sim, A=x[1], B=x[2], d=x[3], bed=bed, max=1e6) } )

# alternativelty, speed improvements can be obtained by optimising A and B for each
# proposed d. Then doing 1-D optimisation of d.
# take mean value of simulated data for starting value of a and guess at B
res2 = do_optimisation(pos, sim, bed, A=mean(sim), B=0.0002)
pred = get_expected(pos, A=res2$par[1],  B=res2$par[2], d=res2$par[3], bed=bed, max=1e6)

# plot simulated data (points) with predictions from inferred parameters (line)
# also plot locations of exons as boxes near x axis.
pdf("example1.pdf", width=10, height=6)
plot(pos, sim, pch=".")
points(pos, pred, col="red", type="l")
apply(bed, 1, addexon, col="red", border="red", ylim=c(0, 0.1))
dev.off()




# Now to fit a model with two different types of element - CNEs and exons.
# here we have two extra parameters (an extra B and d for the second element type).
# read in bed files
exonbed = read.bed("sample/exons.bed")
cnebed = read.bed("sample/cnes.bed")


# simulate some diversity data using our model
pos = seq(2254500, 5089500, 1000)
par = c(6.555099e-02, 2.680702e-05, 8.494293e-05, 3.860168e+04, 7.186175e+03)
sim = get_expected_comb(pos, par, bed1=exonbed, bed2=cnebed, max1=1e5, max2=1e5)
act = sim
sim = sim+rnorm(length(sim), mean=0, sd=0.005)

# do 2D optimisation on d1 & d2 (optimising other paramters at each stage.
res = do_optimisation_comb(pos, sim, exonbed, cnebed, par=par, max=1e5)
pred = get_expected_comb(pos, res$par, bed1=exonbed, bed2=cnebed, max1=1e5, max2=1e5)


pdf("example2.pdf", width=10, height=6)
plot(pos, sim, pch=".")
points(pos, pred, col="red", type="l")
points(pos, act, col="blue", type="l")
apply(cnebed, 1, addexon, col="green", border="green", ylim=c(0.01, 0.015))
apply(exonbed, 1, addexon, col="red", border="red", ylim=c(0.01, 0.015))
dev.off()

