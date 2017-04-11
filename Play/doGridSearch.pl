#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long;

my $self=1;


my $rscript = 'Rscript scripty.R';

my @neutralDiversity = ();
my @reductionPerElement = ();
my @recoveryRate =();
for (my $i=1; $i <= 10; $i+=1) {

    push(@neutralDiversity, $i/100);
    push(@reductionPerElement, 1/10**$i);
    push(@recoveryRate, 10**$i);
}

if ($self) {
    foreach my $neu (@neutralDiversity) {
	foreach my $reduction (@reductionPerElement) {
	    foreach my $recovery (@recoveryRate) {
		print "$rscript $neu $reduction $recovery $neu $reduction $recovery\n";
	    }
	}
    }
} else {
    foreach my $neu (@neutralDiversity) {
	foreach my $reduction (@reductionPerElement) {
	    foreach my $recovery (@recoveryRate) {
#		foreach my $neu2 (@neutralDiversity) {
		    foreach my $reduction2 (@reductionPerElement) {
#			foreach my $recovery2 (@recoveryRate) {
			    print "$rscript $neu $reduction $recovery $neu $reduction2 $recovery\n";
			}
#		    }
#		}
	    }
	}
    }

}





