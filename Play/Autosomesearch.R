source("lib/functions.R")
source("lib/single-chr.R")

# load compiled C routines.
dyn.load("c-routines/do_exp_sums.so")
dyn.load("c-routines/calc_rel_dist.so")

args <- commandArgs(T)

if (is.na(args[5]) ) stop("Gimme some parameters!")
par <- as.numeric (args[1:5])

exonBed = read.bed.gmap("primateElements.inProteinExon")
nonexonBed = read.bed.gmap("primateElements.nonProteinExon")

chroms <- 1:22

exonBedList = list()
nonexonBedList = list()
pidList = list()
posList = list()

x <- 1
for (i in chroms) {

    chrom = paste('chr', i, sep='')

    dataFile <- paste("2kb/", i, ".wga", sep='')
    tbl <- read.table(dataFile)
    tbl <- tbl[tbl$V3 >= 1000,]
    tbl <- tbl[tbl$V19 > 0.0,]
    tbl <- within(tbl, piD <- (V8/V3)/V19)
    q <- quantile(tbl$piD, 0.999)

    tbl$V22 <- tbl$V22 * 1e6
    tbl$V21 <- tbl$V21 * 1e6
    tbl <- tbl[ tbl$piD <= q,]

    exonBedList[[x]] <- exonBed[exonBed$chrom == chrom,]
    nonexonBedList[[x]] <- nonexonBed[nonexonBed$chrom == chrom,]

    tbl <- tbl[tbl$V22 < max(nonexonBedList[[x]]$gEnd),]
    tbl <- tbl[tbl$V21 > min(nonexonBedList[[x]]$gStart),]
    tbl <- tbl[tbl$V22 < max(exonBedList[[x]]$gEnd),]
    tbl <- tbl[tbl$V21 > min(exonBedList[[x]]$gStart),]
    
    pidList[[x]] <- tbl$piD
    posList[[x]] <- ((tbl$V22 + tbl$V21 )/ 2)
    x <- x + 1
}


best2<- optim(par, function(y) { 
s <- 0
x <- 1
for (i in chroms) {
    s <- s + sum_sq_comb_gmap(posList[[x]], pidList[[x]], bed1=exonBedList[[x]], bed2=nonexonBedList[[x]], y, mult=T, max1=1e6,
max2=1e6,trunc=T)
    x <- x + 1
}

s
})


print( c(best2$value, best2$par))

