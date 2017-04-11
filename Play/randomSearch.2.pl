#!/usr/bin/env perl


use strict;
use Math::Random;
use Getopt::Long;

my $chromosome='';
my $autosomes=0;
my $rscriptOptions='';
my $useRegression=0;

GetOptions('c=s' => \$chromosome,
	   'r=s' => \$rscriptOptions,
	   'u' => \$useRegression, # guesses the shape parameter based off of a regression estimate
	   'a' => \$autosomes);

my %bounds = qw(
p1low 0.01
p1high 0.5
p2low 0.0000001
p2high 0.001
p3low 0.000000001
p3high 0.001
p4low 100
p4high 10000000
p5low 100
p5high 10000000
);

my $log10 = log(10);
my @chroms = (1..22);
push(@chroms, 'X');

if ($chromosome) {
    @chroms = ($chromosome);
}

my $command = "Rscript stoopidChromsearch.test2.R $rscriptOptions";
my $autosCommand = "Rscript stoopidAutosomesearch.R $rscriptOptions";

while (1) {
    my @params = ();

    if ($useRegression) {
	push(@params, 
	     random_uniform(1, 0.06, 0.08));
    } else {
	push(@params, 
	     random_uniform(1, $bounds{p1low}, $bounds{p1high}));
    }
# loguniform
    push(@params, 
	 10**random_uniform(1, log($bounds{p2low})/$log10, log($bounds{p2high})/$log10));
    push(@params, 
	 10**random_uniform(1, log($bounds{p3low})/$log10, log($bounds{p3high})/$log10));

    if ($useRegression) {
	push(@params, exp(log($params[1])*-1.96262 - 9.95253)); # parameters found from a linear regression on the log-log transform of the shape (p4) vs selection p2) parameters from the top 1000 best hits on chromosome 22
	push(@params, exp(log($params[2])*-1.65903 - 6.88343)); # parameters found from a linear regression on the log-log transform of the shape (p3) vs selection p5) parameters from the top 1000 best hits on chromosome 22
    } else {
	push(@params, 
	     10**random_uniform(1, log($bounds{p4low})/$log10, log($bounds{p4high})/$log10));
	push(@params, 
	     10**random_uniform(1, log($bounds{p5low})/$log10, log($bounds{p5high})/$log10));
    }

    my @results = ();

    if ($autosomes) {
	my $ans = `$autosCommand @params`;
	$ans =~ s/\s+/ /g;
	print "@params\n";
	print 'Autos '  , $ans , "\n";
    } else {
	foreach (@chroms) {
#	    my $ans = `$command chr$_ 2kb/$_.wga @params`;
#	    die "$command -c chr$_ @params";
	    my $ans = `$command -c chr$_ @params`;
	    $ans =~ s/\s+/ /g;
	    if (length($ans) > 4) {
		print "@params\n";
		print $_ , ' '  , $ans , "\n";
	    }
	}
    }
}



