# The human neutralome

### This is a modified version of Halligan's Model C of: (Halligan et al. 2013 : Contributions of Protein-Coding and Regulatory Change to Adaptive Molecular Evolution in Murid Rodents)

The code contained in this repository is far from general use-- it was constructed to compute an answer to the question at hand-- that is, what is the predicted diversity at a locus conditioned on the local landscape of conserved sequences. Further, much of the code is not general-- it is specific to humans and other nonhuman great apes, and it uses static names of files and directories that are best left unchanged. Further, much of the "meat" of this code was written by Daniel Halligan. I (August Woerner) merely tweaked it, augmenting it to take in more direct information on the genetic map. Functions that were modified by me have been clearly outlined, with commenting to describe the changes that I did.

In principal this code requires:
Two sets of conserved elements with filenames:

    primateElements.inProteinExon
    primateElements.nonProteinExon

Note that the *names* of these files may or may not correspond to what's actually in them. Nevertheless, these names are not mutable. These files must be in bed format, and must use UCSC-style chromosome naming (e.g., the name of chromosome 1 is *chr1* not *1*). As per the bed format specification, the first 3 columns must correspond to a conserved element's location. 2 more columns are expected: the start and end position of the element in the genetic map.

You also must have the *lib* directory, the *c-routines* directory, and the *10kb* directory. The *lib* directory contains the R code used in these scripts. The *c-routines* directory contains the c-code that makes parameter optimization slightly less computationally expensive. You will need to run the command `make` in the *c-routines* directory before using this code. 

The *10kb* directory is the place where the information on the observed diversity and divergence values are kept. Initially I used 10kb loci, though again 10kb is just a name; any locus size is supported. Each chromosome has it's own file; and note the lack of the 'chr' in the filenames. For the data in humans these files are in ".wga" format, though the data in nonhuman great apes I switched this to bed format just to make it more palatable and general use. As the .wga files are a little cryptic, I will instead describe what is required for the .bed versions. For the .bed files you need the chromosome name, start and stop position (as per the bed specification). Next what is required is the nucleotide diversity of the locus, and the number of bases used to estimate diversity (to compute pi/bp). Next is the divergence (total), and the number of bases used to estimate divergence (ideally the same as with diversity), again with the ratio of these two variables being necessary to compute divergence/bp. The last two columns are the positions in the genetic map of the locus (start and stop coordinates in genetic units (cM)).
 
Once these files are in place, I have three R scripts that serve to perform a single iteration of [Nelder-Mead](https://en.wikipedia.org/wiki/Nelder%E2%80%93Mead_method) optimization using Model C of [Halligan et al.](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1003995). In brief, we have our observed diversity (really pi/D), and we seek to find 5 parameters that describe a diversity trough around two types of conserved elements (exonic and nonexonic). We then pick those parameters that minimize the sum of squares (SS) difference between the observed diversity and that predicted by our model.

There are some difficulties in performing this computation in practice. First, Nelder-mead optimization is slow, and we need to run it on the entirety of the genome. Further, Nelder-Mead is used to find local opimima. To accommodate the local optimum problem I simply tested many start locations in this search. To do this I used a very simple randomized strategy; pick 5 points at "random" (see *randomSearch.pl* to see what exactly that means), do Nelder Mead, and take the best result (parameters that yielded the lowest SS error). This worked fine for the X chromosome, but it's only ~150 Mb. For the Autosomes I used a two-pass heuristic; first sample 5% of sites at random (but deterministically), and then pick 5 random starting points. Take the 1,000 or 10,000 best starting parameters, and then use these as starting positions for the same search on the entirety of the autosomes.

To implement a single search of, say, chromosome X type:
   
   perl randomSearch.pl -c X

This will perform a single Nelder Mead optimization and print 
the five parameters from Model C to standard out. In particular, it'll print two lines of output, the first will be the starting coordinates, and the second will be the inferred parameters (ignore the negative signs; the search is based on their absolute values).

For a random search of the autosomes (5% randomly sampled loci) do:

    perl randomSearch.pl -s -a

(-s stands for sample, -a stands for autosomes). Will do the same search-style on the autosomes. From this you'll need to parse out the best SS values, and rerun them using the whole of the autosomes. For an example of how to do that see Rerun.example (this is a simple shell script).

Note that pretty much the entirety of this program relies on having a large computer cluster to distribute the work upon. If that's not present, doing this sort of analysis at a genomic scale isn't recommended. I have helper submit scripts for the portable batch system (PBS) in the Play directory. They perhaps will help with this...

If you have any questions please contact me, August, at
August.Woerner AT unthsc DOT edu

 
 




