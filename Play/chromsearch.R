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
dataFile <- '2kb/22.wga'
if (! is.na(args[1])) chrom <- args[1]
if (! is.na(args[2])) dataFile <- args[2]

bed <- gbed[ gbed$chrom == chrom,]

tbl <- read.table(dataFile)

tbl <- tbl[tbl$V3 >= 1000,]
tbl <- tbl[ tbl$V19 > 0.0,]
tbl <- within(tbl, piD <- (V8/V3)/V19)
q <- quantile(tbl$piD, 0.999)


tbl$V22 <- tbl$V22 * 1e6
tbl$V21 <- tbl$V21 * 1e6

tbl <- tbl[ tbl$piD <= q,]
tbl <- tbl[tbl$V22 < max(bed$gEnd),]
tbl <- tbl[tbl$V21 > min(bed$gStart),]


sim <- tbl$piD
pos <- ((tbl$V22 + tbl$V21 )/ 2)

a <- quantile(sim, 0.95)
b <- quantile(sim, 0.05)


effElemLength<-10000

p3Min <- 2^5
p3Max <- 2^30

stop <- 20

results <- data.frame(matrix(ncol = 4, nrow = stop))
colnames(results) <- c('sumsq', 'p1', 'p2', 'p3')
x <- 1

p1 <- a
p2Start <- log(b/p1)/ -effElemLength


for (i in 0:stop) {
  p1 <- a + a*i/(stop*2)
  p2 <- log(b/p1)/ -(effElemLength*(i+1))
  fp <- function(x) { res = get_expected_gmap(pos, A=p1,  B=p2, d=x, bed=bed, max=1e6, mult=T, trunc=T); sum( (res-sim)^2) }
  best <- optimize(fp, c(p3Min, p3Max),tol=1)
 # best <- do_optimisation_gmap(pos, sim, bed, A=p1, B=p2)
 # par <- best$par
  par <- c(p1,p2, best$minimum) 
  best2<- optim(par, function(x) { sum_sq_gmap(pos, sim, A=x[1], B=x[2],d=x[3], bed=bed, max=1e6,mult=T) } ) 
 
  results[x,] <- c(best2$value, best2$par[1], best2$par[2], best2$par[3])
  print(results[x,])
  
  x<- x + 1   
}


exonBed = read.bed.gmap("primateElements.inProteinExon")
nonexonBed = read.bed.gmap("primateElements.nonProteinExon")
exonBed <- exonBed[ exonBed$chrom == chrom,]
nonexonBed <- nonexonBed[ nonexonBed$chrom == chrom,]

par <- c(result[3,2],  result[3,3], result[3,3], result[3,4], result[3,4] ) 
best2<- optim(par, function(x) { sum_sq_comb_gmap(pos, sim, bed1=exonBed, bed2=nonexonBed, x, mult=T, max1=1e6, max2=1e6)})

#get_expected_comb_gmap = function(pos, par, bed1, bed2, max1, max2, mult=T, bgsel=F, trunc=T, bgpred=NULL, ...) {

fit <- get_expected_comb_gmap(pos, best2$par, exonBed, nonexonBed, max1=1e6, max2=1e6, mult=T, trunc=T)

protexonBed <- read.bed.gmap("primateElementsInExonsUnionExons")
protexonBed <- protexonBed[ protexonBed$chrom == chrom,]
best3<- optim(par, function(x) { sum_sq_comb_gmap(pos, sim, bed1=protexonBed, bed2=nonexonBed, x, mult=T, max1=1e6, max2=1e6)})
print(c(best3$value, best3$par))

#fit3 <- get_expected_comb_gmap(pos, best3$par, protexonBed, nonexonBed, max1=1e6, max2=1e6, mult=T, trunc=T)
