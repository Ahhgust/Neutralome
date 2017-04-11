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
  opt$locusdir='2kb'
if (is.null(opt$chrom))
   opt$chrom = 'chr22'
if (is.null(opt$mult)) 
  opt$mult=F


chrom = opt$chrom;

#dataFile <- '2kb/22.wga'
dataFile <- paste(opt$locusdir , "/", substr(chrom, 4, 7), ".wga", sep='')



exonBed = read.bed.gmap('primateElements.inProteinExon' )#"primateElements.inProteinExon")
nonexonBed = read.bed.gmap("primateElements.nonProteinExon.encodeStrongActive.cM")
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
tbl <- tbl[tbl$V22 < max(nonexonBed$gEnd),]
tbl <- tbl[tbl$V21 > min(nonexonBed$gStart),]
tbl <- tbl[tbl$V22 < max(exonBed$gEnd),]
tbl <- tbl[tbl$V21 > min(exonBed$gStart),]


sim <- tbl$piD
if (opt$physical) {
stop("broken")
# very annoying code to parse chr22:1-200 into 1 and 200
  tbl <- within(tbl, start <- as.numeric(strsplit( as.character(V1),
  split="[:-]",  perl=T)[[1]][2]))
  tbl <- within(tbl, stop <- as.numeric(strsplit( as.character(V1),
  split="[:-]",  perl=T)[[1]][3]))
  pos <- ( (tbl$start + tbl$stop)/2)
  best2<- optim(par, function(x) { sum_sq_comb(pos, sim, bed1=exonBed,
  bed2=nonexonBed, x, mult=opt$mult, max1=1e6, max2=1e6,trunc=T)})
} else {
  pos <- ((tbl$V22 + tbl$V21 )/ 2)# position is in the genetic map
  best2<- optim(par, function(x) { sum_sq_comb_gmap(pos, sim, bed1=exonBed,
  bed2=nonexonBed, x, mult=opt$mult, max1=1e6, max2=1e6,trunc=T)})
}

print( c(best2$value, best2$par))

