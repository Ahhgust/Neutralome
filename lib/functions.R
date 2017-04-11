################################################################################
# function to read (& sort) a bed file
################################################################################
read.bed = function(file, sort=F, strandAsNumeric=F, ...) {
   x = read.table(file, stringsAsFactors=F, ...)
   names = c("chrom","chromStart","chromEnd","name","score","strand",
      "thickStart","thickEnd","itemRgb","blockCount","blockSizes","blockStarts")
   colnames(x) = names[1:ncol(x)]
   if(is.null(x$name))   { x$name="undef"}
   if(is.null(x$score))  { x$score="undef"}
   if(is.null(x$strand)) { x$strand="+"}
   x$chromStart = as.numeric(x$chromStart)
   x$chromEnd = as.numeric(x$chromEnd)
   if(sort) { 
     x = x[order(x$chrom, x$chromStart),]
   }
   if(strandAsNumeric) {
      x$strand[which(x$strand == '+')] = 1
      x$strand[which(x$strand == '-')] = -1
      x$strand = as.numeric(x$strand)
   }
   x
}

# Written by August Woerner
# This reads in a bed file format as read.bed
# except that this bed file has been augmented to also included positions in
# the genetic maps
read.bed.gmap = function(file, sort=F, strandAsNumeric=F, ...) {
   x = read.table(file, stringsAsFactors=F, ...)
   names = c("chrom","chromStart","chromEnd","gStart","gEnd","strand",
      "thickStart","thickEnd","itemRgb","blockCount","blockSizes","blockStarts")
   colnames(x) = names[1:ncol(x)]
   if(is.null(x$name))   { x$name="undef"}
   if(is.null(x$score))  { x$score="undef"}
   if(is.null(x$strand)) { x$strand="+"}
   

   x$gStart = as.double(x$gStart)
   x$gEnd = as.double(x$gEnd)
   x$gStart <- x$gStart * 1000000.0; # scale the genetic map such that the genetic units
   x$gEnd <- x$gEnd * 1000000.0; # are on approximately the same scale as the physical units

   x$chromStart = as.numeric(x$chromStart)
   x$chromEnd = as.numeric(x$chromEnd)
   if(sort) { 
     x = x[order(x$chrom, x$chromStart),]
   }
   if(strandAsNumeric) {
      x$strand[which(x$strand == '+')] = 1
      x$strand[which(x$strand == '-')] = -1
      x$strand = as.numeric(x$strand)
   }
   x
}



##########################################################################################
# Pure R function to calculate predicted relative reduction in diversity 
# at a single position given a data frame (ann) containing start and end positions
# of neighbouring conserved blocks.
# d is a distance parameter (from an expoential distribution)
# max = maximum distance to consider around position
##########################################################################################
do_sum = function(pos, ann, d, max=1e5) {
   d1 = abs(ann$start - pos)
   d2 = abs(ann$end - pos)
   ann_sub = ann[which(d1 <= max | d2 <= max),]
   if(nrow(ann_sub)==0) { return(0) }
   ann_exp = unlist(apply(ann_sub, 1, function(r) { seq((r[1]+1), r[2], 1) }))
   d3 = abs(ann_exp - pos)
   sum(exp(-d3/d))
}

##########################################################################################
# Uses call to a c routine to calculate predicted relative reduction in diversity 
# at a single position given a data frame (ann) containing start and end positions
# Replicates "do_sum".
##########################################################################################
do_sum_c = function(pos, ann, d, max=1e5) {
   d1 = abs(ann$start - pos)
   d2 = abs(ann$end - pos)
   ann_sub = ann[which(d1 <= max | d2 <= max),]
   if(nrow(ann_sub)==0) { return(0) }
   .C("do_exp_sum", 
   		as.integer(pos), 
   		as.integer(nrow(ann_sub)), 
   		as.integer(ann_sub$start), 
   		as.integer(ann_sub$end), 
   		as.double(d),
   		as.double(0.0)
   	)[[6]]
}

# calculate sum of negative exponential expectations.
# if ignore == T, data points with no features within max distance are given 
# value NA rather than the value of A
do_exp_sums = function(pos, bed, d, max=1e6, weighted=F, ignore=F, s='chromStart', e='chromEnd') {
   out <- .C("do_exp_sums", 
      as.integer(pos),
      as.integer(length(pos)), 
      as.integer(nrow(bed)), 
      as.integer(bed[[s]]), 
      as.integer(bed[[e]]), 
      as.double(d), 
      as.double(max), 
      as.integer(weighted), 
      as.integer(ignore), 
      result = double(length(pos))
   )
   i = which(out[['result']]==-1)
   out[['result']][i] = NA
   out[['result']]
}

# calculate sum of negative exponential expectations.
# if ignore == T, data points with no features within max distance are given 
# value NA rather than the value of A
# Written by August Woerner
# Again, analog to do_exp_sums but instead uses positions in the genetic maps
do_exp_sums_gmap = function(pos, bed, d, max=1e6, weighted=F, ignore=F, s='chromStart', e='chromEnd',
gs='gStart', ge='gEnd'
) {
   out <- .C("do_exp_sums_gmap", 
      as.double(pos),
      as.integer(length(pos)), 
      as.integer(nrow(bed)), 
      as.integer(bed[[s]]), 
      as.integer(bed[[e]]), 
      as.double(bed[[gs]]), 
      as.double(bed[[ge]]), 
      as.double(d), 
      as.double(max), 
      as.integer(weighted), 
      as.integer(ignore), 
      result = double(length(pos))
   )
   i = which(out[['result']]==-1)
   out[['result']][i] = NA
   out[['result']]
}

# calculate distance to nearest element border from mid point of non-overlapping window
calc_rel_dist = function(data, cbed, lwr, upr, p='pos', s='chromStart', e='chromEnd') {
	lwr[which(is.na(lwr))] = -1
	upr[which(is.na(upr))] = -1
    out <- .C("calc_rel_dist", 
      as.double(data[[p]]),
      as.integer(nrow(data)), 
      as.double(lwr), 
      as.double(cbed[[s]]), 
      as.double(cbed[[e]]), 
      as.double(upr), 
      as.integer(cbed$strand), 
      as.integer(nrow(cbed)), 
      dist = double(nrow(data)),
      len = double(nrow(data))
   )
   dist = out[['dist']]
   len = out[['len']]
   dist[which(dist==0)] = NA
   len[which(len==0)] = NA
   data.frame(dist = dist,len = len)
}


# calculate sum of negative exponential expectations.
# if ignore == T, data points with no features within max distance are given 
# value NA rather than the value of A
do_exp_sums = function(pos, bed, d, max=1e5, weighted=F, ignore=F, s='chromStart', e='chromEnd') {
   out <- .C("do_exp_sums", 
      as.integer(pos),
      as.integer(length(pos)), 
      as.integer(nrow(bed)), 
      as.integer(bed[[s]]), 
      as.integer(bed[[e]]), 
      as.double(d), 
      as.double(max), 
      as.integer(weighted), 
      as.integer(ignore), 
      result = double(length(pos))
   )
   i = which(out[['result']]==-1)
   out[['result']][i] = NA
   out[['result']]
}



exp_bg_red = function(pos, bed, vec, max=1e5, s='chromStart', e='chromEnd') {
   out <- .C("exp_bg_red", 
      as.double(pos),
      as.double(length(pos)), 
      as.double(nrow(bed)), 
      as.double(bed[[s]]), 
      as.double(bed[[e]]), 
      as.double(vec), 
      as.double(length(vec)), 
      as.double(max), 
      result = double(length(pos))
   )
   out[['result']]
}

# multiplicative model of above
exp_bg_red_m = function(pos, bed, vec, max=1e5, s='chromStart', e='chromEnd') {
   out <- .C("exp_bg_red_m", 
      as.double(pos),
      as.double(length(pos)), 
      as.double(nrow(bed)), 
      as.double(bed[[s]]), 
      as.double(bed[[e]]), 
      as.double(vec), 
      as.double(length(vec)), 
      as.double(max), 
      result = double(length(pos))
   )
   out[['result']]
}

