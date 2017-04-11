
# mult=T calculates an approximation to a multiplicative model.
# trunc = T means that negative exp values are truncated such that minimum==0
get_expected = function(pos, A, B, d, trunc=T, mult=T, ...) {
   exp = vector()
   if(mult) {
      exp = exp(log(A) - B*do_exp_sums(pos, d=d, ...))
   } else {
      exp = A * (1 - B*do_exp_sums(pos, d=d, ...))
   }      
   if(trunc) {
     exp[which(exp<0)] = 0
   }
   exp   
}

# mult=T calculates an approximation to a multiplicative model.
# trunc = T means that negative exp values are truncated such that minimum==0
get_expected_gmap = function(pos, A, B, d, trunc=F, mult=T, ...) {

   exp = vector()
   if(mult) {
      exp = exp(log(A) - B*do_exp_sums_gmap(pos, d=d, ...))
   } else {
      exp = A * (1 - B*do_exp_sums_gmap(pos, d=d, ...))
   }      
   
   if(trunc) {
     exp[which(exp<0)] = 0
   }
   exp   
}

sum_sq = function(pos, obs, A, B, d, ...) {
   exp = get_expected(pos, A, B, d, ...)
   sumsq = sum((exp-obs)^2, na.rm=T)
 #  print(c(A, B, d, sumsq))
   sumsq
}

sum_sq_gmap = function(pos, obs, A, B, d, ...) {
   sumsq <- 2/0.0
   
   if (A > 0) {
     if (B>0) {
       if (d > 0) {
       
    	 exp = get_expected_gmap(pos, A, B, d, ...)
      	 sumsq = sum((exp-obs)^2, na.rm=T)
       }
     }
   }

   sumsq
 
}

# function to return sumsq, optimised for A and B
optimise_A_B = function(obs, expsum, par, mult=F) {
   optim(par, function(x) { 
      exp = vector()
      if(mult) {
         exp = exp(log(x[1]) - x[2]*expsum)
      } else {
         exp = x[1] * (1 - x[2]*expsum)
      }      
      sum((exp-obs)^2, na.rm=T)
   })
}

# 1-D optimisation on d
# lower and upper specify the range over which d is optimised.
do_optimisation = function(pos, obs, bed, A, B, lower=0, upper=5e6, mult=T, ...) {
   par = c(A,B)
   d = optimise(function(d) { 
      res = optimise_A_B(obs, do_exp_sums(pos, d=d, bed=bed, ...), par, mult=mult)
#      print(c(res$par, d, res$value))
      res$value
   }, lower=lower, upper=upper)$minimum
   res = optimise_A_B(obs, do_exp_sums(pos, d=d, bed=bed), par, mult=mult)
   list(par=c(res$par, d), value=res$value)
}

# 1-D optimisation on d
# lower and upper specify the range over which d is optimised.
do_optimisation_gmap = function(pos, obs, bed, A, B, lower=0, upper=5e6, mult=T, ...) {
   par = c(A,B)
   d = optimise(function(d) { 
      res = optimise_A_B(obs, do_exp_sums_gmap(pos, d=d, bed=bed, ...), par, mult=mult)
#      print(c(res$par, d, res$value))
      res$value
   }, lower=lower, upper=upper)$minimum
   res = optimise_A_B(obs, do_exp_sums(pos, d=d, bed=bed), par, mult=mult)
   list(par=c(res$par, d), value=res$value)
}


get_expected_comb_gmap = function(pos, par, bed1, bed2, max1, max2, mult=T, bgsel=F, trunc=T, bgpred=NULL, ...) {
   par <- abs(par)
   expsum1 = do_exp_sums_gmap(pos, d=par[4], bed=bed1, max=max1, ...)
   expsum2 = do_exp_sums_gmap(pos, d=par[5], bed=bed2, max=max2, ...)
	if(mult) {
		exp = exp(log(par[1]) - par[2]*expsum1 - par[3]*expsum2)
	} else {
		exp = par[1] * (1 - par[2]*expsum1 - par[3]*expsum2)
	}
	if(bgsel) {
		exp = exp * bgpred
	}
	if(trunc) {
		exp[which(exp<0)] = 0
	}
	return(exp)    
}

sum_sq_comb_gmap = function(pos, obs, bed1, bed2, par, ...) {
   exp = get_expected_comb_gmap(pos, par, bed1, bed2, ...)
   sum((exp-obs)^2, na.rm=T)
}

get_expected_comb = function(pos, par, bed1, bed2, max1, max2, mult=F,
bgsel=F, trunc=F, bgpred=NULL, ...) {
   par <- abs(par)	       
   expsum1 = do_exp_sums(pos, d=par[4], bed=bed1, max=max1, ...)
   expsum2 = do_exp_sums(pos, d=par[5], bed=bed2, max=max2, ...)
	if(mult) {
		exp = exp(log(par[1]) - par[2]*expsum1 - par[3]*expsum2)
	} else {
		exp = par[1] * (1 - par[2]*expsum1 - par[3]*expsum2)
	}
	if(bgsel) {
		exp = exp * bgpred
	}
	if(trunc) {
		exp[which(exp<0)] = 0
	}
	return(exp)    
}

sum_sq_comb = function(pos, obs, bed1, bed2, par, ...) {
   exp = get_expected_comb(pos, par, bed1, bed2, ...)
   sum((exp-obs)^2, na.rm=T)
}


# function to return sumsq, optimised for A1, B1 and B2 under the combined model
# parameters must be in order A, B1, B2
optimise_A_B1_B2 = function(obs, par, expsum1, expsum2, mult=F) {
   optim(par, function(x) { 
      exp = vector()
      if(mult) {
         exp = exp(log(x[1]) - x[2]*expsum1 - x[3]*expsum2)
      } else {
         exp = x[1] * (1 - x[2]*expsum1 - x[3]*expsum2)
      }      
      sum((exp-obs)^2, na.rm=T)
   })
}

optimise_A_B1_B2_bgsel = function(obs, par, expsum1, expsum2, bgred1, bgred2, mult=F) {
   optim(par, function(x) { 
   	exp = vector()
   	if(mult) {
         exp = exp(log(x[1]) - x[2]*expsum1 - x[3]*expsum2) * bgred1 * bgred2
   	} else {
         exp = x[1] * (1 - x[2]*expsum1 - x[3]*expsum2) * bgred1 * bgred2
         exp[which(exp<0)] = 0
      }
      sum((exp-obs)^2, na.rm=T)
   })
}

# 2-D optimisation on d1 and d2 under a combined model
# parameters must be in order A, B1, B2, d1, d2
do_optimisation_comb = function(pos, obs, bed1, bed2, par, mult=F, ...) {
   d_res = optim(par[4:5], function(d) {
      expsum1 = do_exp_sums(pos, d=d[1], bed=bed1, ...)
      expsum2 = do_exp_sums(pos, d=d[2], bed=bed2, ...)
      ab_res = optimise_A_B1_B2(obs, par[1:3], expsum1, expsum2, mult=mult)
      print(c(ab_res$par, d[1], d[2], ab_res$value))
      ab_res$value
   })
	# do final optimisation of A, B1, B2
	expsum1 = do_exp_sums(pos, d=d_res$par[1], bed=bed1, ...)
	expsum2 = do_exp_sums(pos, d=d_res$par[2], bed=bed2, ...)
   ab_res = optimise_A_B1_B2(obs, par[1:3], expsum1, expsum2, mult=mult)
   list(par=c(ab_res$par,d_res$par), value=ab_res$value)
}


# 2-D optimisation on d1 and d2 under a combined model with bgsel
# parameters must be in order A, B1, B2, d1, d2
do_optimisation_comb_bgsel = function(pos, obs, bed1, bed2, par, bgred1, bgred2, ...) {
   d_res = optim(par[4:5], function(d) {
      expsum1 = do_exp_sums(pos, d=d[1], bed=bed1, ...)
      expsum2 = do_exp_sums(pos, d=d[2], bed=bed2, ...)
      ab_res = optimise_A_B1_B2_bgsel(obs, par[1:3], expsum1, expsum2, bgred1, bgred2)
      print(c(ab_res$par, d[1], d[2], ab_res$value))
      ab_res$value
   })
	# do final optimisation of A, B1, B2
	expsum1 = do_exp_sums(pos, d=d_res$par[1], bed=bed1, ...)
	expsum2 = do_exp_sums(pos, d=d_res$par[2], bed=bed2, ...)
   ab_res = optimise_A_B1_B2_bgsel(obs, par[1:3], expsum1, expsum2, bgred1, bgred2)
   list(par=c(ab_res$par,d_res$par), value=ab_res$value)
}
