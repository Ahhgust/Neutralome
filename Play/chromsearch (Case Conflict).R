source("lib/functions.R")
source("lib/single-chr.R")

# load compiled C routines.
dyn.load("c-routines/do_exp_sums.so")
dyn.load("c-routines/calc_rel_dist.so")



chrom = 'chrX';
args <- commandArgs(T)
dataFile <- '10kb/X.wga'
if (! is.na(args[1])) 
if (! is.na(args[2])) 

if (is.na(args[7]) ) stop("Gimme some parameters!")
chrom <- args[1]
dataFile <- args[2]
par <- as.numeric (args[3:7])

exonBed = read.bed.gmap("primateElements.inProteinExon")
nonexonBed = read.bed.gmap("primateElements.nonProteinExon")
exonBed <- exonBed[ exonBed$chrom == chrom,]
nonexonBed <- nonexonBed[ nonexonBed$chrom == chrom,]


tbl <- read.table(dataFile)

tbl <- tbl[tbl$V3 >= 1000,]
tbl <- tbl[ tbl$V19 > 0.0,]
tbl <- within(tbl, piD <- (V8/V3)/V19)
q <- quantile(tbl$piD, 0.999)


tbl$V22 <- tbl$V22 * 1e6
tbl$V21 <- tbl$V21 * 1e6

tbl <- tbl[ tbl$piD <= q,]
tbl <- tbl[tbl$V22 < max(nonexonBed$gEnd) - 1e6,]
tbl <- tbl[tbl$V21 > min(nonexonBed$gStart)+1e6,]
tbl <- tbl[tbl$V22 < max(exonBed$gEnd)-1e6,]
tbl <- tbl[tbl$V21 > min(exonBed$gStart)+1e6,]


sim <- tbl$piD
pos <- ((tbl$V22 + tbl$V21 )/ 2)

best2<- optim(par, function(x) { sum_sq_comb_gmap(pos, sim, bed1=exonBed, 
bed2=nonexonBed, x, mult=T, max1=1e6, max2=1e6,trunc=T)},
control=list(maxit=10000) )



print( c(best2$value, best2$par, best2$convergence))


