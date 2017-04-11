# compile .c code to .so
all:
	R CMD SHLIB do_exp_sums.c
	R CMD SHLIB calc_rel_dist.c

clean:
	rm calc_rel_dist.o  calc_rel_dist.so  do_exp_sums.o  do_exp_sums.so
