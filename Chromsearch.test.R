source("lib/functions.R")
source("lib/single-chr.R")

# load compiled C routines.
dyn.load("c-routines/do_exp_sums.so")
dyn.load("c-routines/calc_rel_dist.so")

library('getopt')
spec <- matrix(c(
'physical', 'p', 0, "logical", # physical or genetic units?
'locusdir', 'd', 1, "character", # the directory for the loci
'chrom', 'c', 1, "character", 
'mult', 'm', 0, "logical" # pseudomultiplicative model
), byrow=T, ncol=4
)

args <- commandArgs(T)
numargs <- length(args)
if (numargs < 5) stop("I need at least 5 arguments!")
par <- as.numeric (args[(numargs-4):numargs])
args <- args[1:(numargs-5)]

opt <- getopt(spec, opt=args)


if (is.null(opt$physical)) 
  opt$physical=F
if (is.null(opt$locusdir))
  opt$locusdir='10kb'
if (is.null(opt$chrom))
   opt$chrom = 'chrX'
if (is.null(opt$mult)) 
  opt$mult=F


chrom = opt$chrom;


dataFile <- paste(opt$locusdir , "/", substr(chrom, 4, 7), ".bed", sep='')



exonBed = read.bed.gmap("primateElements.inProteinExon")
nonexonBed = read.bed.gmap("primateElements.nonProteinExon")
exonBed <- exonBed[ exonBed$chrom == chrom,]
nonexonBed <- nonexonBed[ nonexonBed$chrom == chrom,]


tbl <- read.table(dataFile)
    colnames(tbl) <- c('chrom', 'start', 'stop', 'piOfLocus', 'piLen', 'divOfLocus', 'dLen', 'cmStart', 'cmStop')
    tbl <- tbl[tbl$piLen >= 1000,]
    tbl <- tbl[tbl$dLen >= 1000,]
    tbl <- tbl[tbl$divOfLocus > 0.0,]
    tbl <- within(tbl, piD <- (piOfLocus/piLen)/(divOfLocus/dLen) )
    q <- quantile(tbl$piD, 0.999)

    tbl$cmStart <- tbl$cmStart * 1e6
    tbl$cmStop <- tbl$cmStop * 1e6
    tbl <- tbl[ tbl$piD <= q,]


    tbl <- tbl[tbl$cmStop < max(nonexonBed$gEnd) - 1e6,]
    tbl <- tbl[tbl$cmStart > min(nonexonBed$gStart) + 1e6,]
    tbl <- tbl[tbl$cmStop < max(exonBed$gEnd)-1e6,]
    tbl <- tbl[tbl$cmStart > min(exonBed$gStart)+1e6,]

sim <- tbl$piD

if (opt$physical) {
# very annoying code to parse chr22:1-200 into 1 and 200
stop("this is wrong")
  tbl <- within(tbl, start <- as.numeric(strsplit( as.character(V1),
  split="[:-]",  perl=T)[[1]][2]))
  tbl <- within(tbl, stop <- as.numeric(strsplit( as.character(V1),
  split="[:-]",  perl=T)[[1]][3]))
  pos <- ( (tbl$start + tbl$stop)/2)
  best2<- optim(par, function(x) { sum_sq_comb(pos, sim, bed1=exonBed,
  bed2=nonexonBed, x, mult=opt$mult, max1=1e6, max2=1e6,trunc=T)})
} else {
  pos <- ((tbl$cmStop + tbl$cmStart)/ 2)# position is in the genetic map
  
  best2<- optim(par, function(x) { sum_sq_comb_gmap(pos, sim, bed1=exonBed,
  bed2=nonexonBed, x, mult=opt$mult, max1=1e6, max2=1e6,trunc=T)},
  control=list(maxit=10000))	     
}


print( c(best2$value, best2$par, best2$convergence))


