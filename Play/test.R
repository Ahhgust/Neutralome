source("lib/functions.R")
dyn.load("c-routines/do_exp_sums.so")
dyn.load("c-routines/calc_rel_dist.so")


# read bed file containing locations of conserved regions.
bed = read.bed("sample/cnes.bed")

# create "annotation": a data frame with just the start and end positions
ann = data.frame(
   start = bed$chromStart,
   end   = bed$chromEnd
)

pos = 5080000
d = 150000

# Pure R calculation for a single position
# First get vector of all conserved positions - for each conserved element, get vector of all positions
ann_exp = unlist(apply(ann, 1, function(r) seq(r[1]+1, r[2], 1) ))

# get distance from position to each conserved base.
dist = abs(ann_exp - pos)

# sum over sites of exp(-dist/d)
sum(exp(-dist/150000))

# Use pure R function to do the same.
do_sum(pos, ann, d, max=1e10) 

# Now use a C function
.C("do_exp_sum", as.integer(pos), as.integer(nrow(ann)), as.integer(ann$start), as.integer(ann$end), as.double(d), as.double(0.0))[[6]]


# Set up a vector of positions.
pos <- seq(min(ann), max(ann), 1000)

# apply R function over positions = slow
# uses approximation to only conisder distances < 1e5 from each position
res = sapply(pos, do_sum, ann=ann, d=15000, max=1e5)

# apply c function over positions = better
res2 = sapply(pos, do_sum_c, ann=ann, d=15000, max=1e5)

# use a C function for whole loop = best
res3 = do_exp_sums(pos, bed, d=15000, max=1e5)


# plot results to show they are more or less the same
plot(res, type="l")
lines(res2, type="l", col="red")
lines(res3, type="l", col="blue")






