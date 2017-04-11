#include <R.h>
#include <stdio.h>
#include <stdlib.h>

// assumes bed and dat are sorted by position
// pos      = midpoint in non-overlapping window
// ldat     = length of pos vector
// lwr      = vector of lower bounds for bed flanks (if pos is between lwr and start this is your closest element)
// start    = vector of bed element start
// end      = vector of bed element ends
// upr      = vector of upper bounds for bed flanks (if pos is between end and upr this is your closest element)
// strand   = vector of bed element strands
// lbed     = length of bed vectors
// rel_dist = return vector of distances to nearest bed element
// rel_len  = return vector of length of nearest bed element
void calc_rel_dist(double *pos, int *ldat, double *lwr, double *start, 
    double *end, double *upr, int *strand, int *lbed, double *rel_dist, double *rel_len)
{
	int di = 0;  // counter for data lines

   // run through bed lines, consider valid lwr and upr sections only
   // lwr and upr bed sections will be non-overlapping.
   for(int bi=0; bi<*lbed; bi++) {
   
		// consider lwr block: lwr..start
		if(lwr[bi] != -1) {
			while(pos[di] < start[bi]) { 
			   if(pos[di] < lwr[bi]) { 
			   	    di++;
					if(di >= *ldat) { return; }
			   	    continue; 
				}
				rel_dist[di] = (pos[di] - start[bi]) * strand[bi];
				rel_len[di] =  (end[bi] - start[bi]);
//				Rprintf("lwr: %d %d    %1.1f %1.1f %1.1f %d %1.1f \n", bi, di, pos[di], lwr[bi], start[bi], strand[bi], rel_dist[di]);
				di++; 
				if(di >= *ldat) { return; }
			}  
		}
		
		// consider upr block: end..upr
		if(upr[bi] != -1) {
			while(pos[di] < upr[bi]) { 
			   if(pos[di] < end[bi]) { 
			   	    di++;
					if(di >= *ldat) { return; }
			   	    continue; 
			   }
				rel_dist[di] = (pos[di] - end[bi]) * strand[bi];
				rel_len[di] =  (end[bi] - start[bi]);
//				Rprintf("upr: %d %d     %1.1f %1.1f %1.1f %d %1.1f \n", bi, di, pos[di], lwr[bi], start[bi], strand[bi], rel_dist[di]);
				di++; 
				if(di >= *ldat) { return; }
			}  
		}
   }
}
