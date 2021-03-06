#include <R.h>
#include <stdio.h>
#include <stdlib.h>

#define MAX 1000000

void do_exp_sum(int *pos, int *n, int *start, int *end, double *d, double *sum)
{
	int i,j;
	*sum = 0.0;
	for(i=0; i<*n; i++) {
	   for(j=start[i]+1; j<=end[i]; j++) {
	      double e = (double)abs(j-*pos);
	      *sum += exp(-e / *d);
		}
	} 
}

// p_ind = index in position vector (pos), with length l
// b_ind = index in bed vectors (start, end), with length n
// c_ind = running index in bed vectors
// "weighted" can be true or false. If true, then contributions from each base
// are weighted by 1/l where l is the length of the element. Therefore, 
// contributions from different exons have equal weights, but contributions from
// different bases within those exons will have different weights...
void do_exp_sums(int *pos, int *l, int *n, int *start, int *end, double *d, 
   double *max, int *weighted, int *ignore, double *sum)
{
	int p_ind,j;
	int c_ind = 0;
	int b_ind = 0;
	int count = 0;  // count of number of selected bases near to a position
	
	// compute vector of negative exponential values for specified d value
	if(*max > MAX) { 
	   Rprintf("Memory error!\n");
	   return;
	}
	double denom = 1 - exp(1.0/ *d);

	double exp_vec[MAX+1];
	for(int i=0; i<=MAX; i++) {
	   exp_vec[i] = exp(((double)-i) / *d);
	}
	*max = MAX; // by AW. Let's just eliminate this variable for now...
   // for each position
	for(p_ind=0; p_ind<*l; p_ind++) {
		// increment index in bed vectors if current position is too low
		if(end[c_ind] < pos[p_ind]) {
		    while((pos[p_ind] - end[c_ind]) > *max) { 
			c_ind++;
		    }
		}
		
		sum[p_ind] = 0.0;
		count = 0;
      
		// for bed each element (from c_ind till bed is out of range)

		for(b_ind=c_ind; b_ind<*n; b_ind++) {
		   if(((double)start[b_ind] - (double)pos[p_ind]) > *max) { 
		   	break; 
		   }
//			Rprintf("%d %d %d %f\n", b_ind, start[b_ind], pos[p_ind], *max);
			// for each base in this element (start is 0-based, end is 1-based)
		   int len = end[b_ind]-start[b_ind];
/*
		   int dist = abs(start[b_ind]+1-pos[p_ind]);
		   if(dist > *max) 
      		       continue; 
*/
/*

		   for(j=start[b_ind]+1; j<=end[b_ind]; j++) {
//				double e = (double)abs(j-pos[p_ind]);
//				Rprintf("%lf\n", e);
//				sum[p_ind] += exp(-e / *d);

		       int dist = abs(j-pos[p_ind]);
		       if(dist > *max) { continue; }
		       if(dist == 0)  { dist = 1; } // hack to avoid not reducing diversity if we are in a CNE / exon
		       count++;
//				Rprintf("%i\n", e);
		       if(*weighted==1) {
			   sum[p_ind] += exp_vec[dist] / len;
		       } else {
			   sum[p_ind] += exp_vec[dist];
//			   sum[p_ind] += exp( ( (double)-dist) / *d); // prohibitively slow
		       }
		   }
*/

		   int dist1 = start[b_ind]+1 - pos[p_ind];
		   int dist2 = end[b_ind] - pos[p_ind];
		   int j,n;
		   if (dist1 < 0 && dist2 < 0) {
		       j = -dist2;
		       n = -dist1;
		   } else if (dist1 > 0 && dist2 > 0) {
		       j = dist1;
		       n = dist2;
		   } else {
// sum over the left-hand side
		       sum[p_ind] += - (( 1 - 
					  exp( dist1/ *d) ) / denom);
		       j = 1;
		       n = dist2;
		   }

// can be done using LUTs as well...
		   sum[p_ind] += -(( exp( (1 - j)/ *d) - 
				    exp( -n/ *d) ) / denom);

		}
      
      if(*ignore==1 && count==0) {
         sum[p_ind] = -1;
      }
      
//		Rprintf("%f %i\n", sum[p_ind], count);
   }
}




// MULTIPLICATIVE MODEL - note that we can't just calculate the sum of 
// exponential values, but must calculate the product of [1 - B.exp(-x/d)]
// (the expected proportional reductions in diversity)

// p_ind = index in position vector (pos), with length l
// b_ind = index in bed vectors (start, end), with length n
// "weighted" can be true or false. If true, then contributions from each base
// are weighted by 1/l where l is the length of the element. Therefore, 
// contributions from different exons have equal weights, but contributions from
// different bases within those exons will have different weights...
void get_exp_red_mult(int *pos, int *l, int *n, int *start, int *end, double *B, 
   double *d, double *max, int *weighted, int *ignore, double *red)
{
	int p_ind,j;
	int c_ind = 0;
	int b_ind = 0;
	int count = 0;  // count of number of selected bases near to a position
	
	// compute vector of negative exponential values using the d parameter provided
	if(*max > MAX) { 
	   Rprintf("Memory error!\n");
	   return;
	}
	double exp_vec[MAX+1];
	for(int i=0; i<=MAX; i++) {
	   exp_vec[i] = exp(((double)-i) / *d);
	}

   // for each position
	for(p_ind=0; p_ind<*l; p_ind++) {
		// increment index in bed vectors if current position is out of range
		if(end[c_ind] < pos[p_ind]) {
			while((pos[p_ind] - end[c_ind]) > *max) { 
				c_ind++;
			}
		}
		
		red[p_ind] = 1.0;
      count = 0;
      
		// for bed each element (from c_ind till bed is out of range)
		for(b_ind=c_ind; b_ind<*n; b_ind++) {
		   if(((double)start[b_ind] - (double)pos[p_ind]) > *max) { 
		   	break; 
			}
//			Rprintf("%d %d %d %f\n", b_ind, start[b_ind], pos[p_ind], *max);
			// for each base in this element (start is 0-based, end is 1-based)
			int len = end[b_ind]-start[b_ind];
			for(j=start[b_ind]+1; j<=end[b_ind]; j++) {
//				double e = (double)abs(j-pos[p_ind]);
//				Rprintf("%lf\n", e);
//				red[p_ind] += exp(-e / *d);

				int dist = abs(j-pos[p_ind]);
			   if(dist > *max) { continue; }
				if(dist == 0)   { continue; }
				count++;
//				Rprintf("%i\n", e);				
//            Rprintf("%lf  %lf  %i  %lf  %lf   %lf\n", *d, *B, dist, exp_vec[dist], (1 - *B * exp_vec[dist]), red[p_ind]);
				if(*weighted==1) {
					red[p_ind] *= 1.0 - *B * (exp_vec[dist] / len);
				} else {
					red[p_ind] *= 1.0 - *B * (exp_vec[dist]);
				}
			}
		}
      
      if(*ignore==1 && count==0) {
         red[p_ind] = -1;
      }
      
//		Rprintf("%f %i\n", red[p_ind], count);
   }
}


//   out <- .C("exp_bg_red", 
//      as.integer(pos),
//      as.integer(length(pos)), 
//      as.integer(nrow(bed)), 
//      as.integer(bed$chromStart), 
//      as.integer(bed$chromEnd), 
//      as.double(vec), 
//      as.double(length(vec)), 
//      as.double(max), 
//      result = double(length(pos))
//   )
void exp_bg_red(double *pos, double *l, double *n, double *start, double *end, 
   double *bgvec, double *veclen, double *max, double *sum)
{
	long p_ind, j, b_ind;
	long c_ind = 0;
	
	for(p_ind=0; p_ind<(long) *l; p_ind++) {
		// increment index in bed vectors if current position is out of range
		if(end[c_ind] < pos[p_ind]) {
			while((pos[p_ind] - end[c_ind]) > *max) { 
				c_ind++;
			}
		}
		
		sum[p_ind] = 0.0;

		// for bed each element (from c_ind till bed is out of range
		for(b_ind=c_ind; b_ind<(long) *n; b_ind++) {
		   if((start[b_ind] - pos[p_ind]) > *max) { 
		   	break; 
			}
//			Rprintf("%d %d %d %f\n", b_ind, start[b_ind], pos[p_ind], *max);
			// for each base in this element (start is 0-based, end is 1-based)
			for(j=(long)start[b_ind]+1; j<=(long) end[b_ind]; j++) {
				long dist = abs(j-(long)pos[p_ind]);
			   if(dist > *max) { continue; }
				if(dist == 0)   { continue; }
				if(dist > (long)*veclen) { 
				   dist = (long)*veclen - 1;
//				   Rprintf("%d %f\n", dist, bgvec[dist]);
            }
				sum[p_ind] += bgvec[dist];
//				Rprintf("%i\n", dist);
				
			}
//			Rprintf("%f %f %f %i\n", start[b_ind], end[b_ind], sum[p_ind], p_ind);
		}
//		Rprintf("%d %f\n", p_ind, sum[p_ind]);
   }
}

// this is a multiplicative version of the above.
void exp_bg_red_m(double *pos, double *l, double *n, double *start, double *end, 
   double *bgvec, double *veclen, double *max, double *prod)
{
	long p_ind, j, b_ind;
	long c_ind = 0;
	
	for(p_ind=0; p_ind<(long) *l; p_ind++) {
		// increment index in bed vectors if current position is out of range
		if(end[c_ind] < pos[p_ind]) {
			while((pos[p_ind] - end[c_ind]) > *max) { 
				c_ind++;
			}
		}
		
		prod[p_ind] = 0.0;

		// for bed each element (from c_ind till bed is out of range
		for(b_ind=c_ind; b_ind<(long) *n; b_ind++) {
		   if((start[b_ind] - pos[p_ind]) > *max) { 
		   	break; 
			}
//			Rprintf("%d %d %d %f\n", b_ind, start[b_ind], pos[p_ind], *max);
			// for each base in this element (start is 0-based, end is 1-based)
			for(j=(long)start[b_ind]+1; j<=(long) end[b_ind]; j++) {
				long dist = abs(j-(long)pos[p_ind]);
			   if(dist > *max) { continue; }
				if(dist == 0)   { continue; }
				if(dist > (long)*veclen) { 
				   dist = (long)*veclen - 1;
//				   Rprintf("%d %f\n", dist, bgvec[dist]);
            }
				prod[p_ind] += log(1 - bgvec[dist]);
//				Rprintf("%d %f\n", dist, 1-bgvec[dist]);
				
			}
//		Rprintf("%f %f %f %i\n", start[b_ind], end[b_ind], prod[p_ind], p_ind);
		}
		prod[p_ind] = exp(prod[p_ind]);
//		Rprintf("%d %f\n", p_ind, sum[p_ind]);
   }
}

