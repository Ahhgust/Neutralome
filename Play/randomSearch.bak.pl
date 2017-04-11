#!/usr/bin/env perl


use strict;
use Math::Random;
use Getopt::Long;

my $chromosome='';
my $autosomes=0;
my $rscriptOptions='';
my $parameterFile = '';

GetOptions('c=s' => \$chromosome,
	   'r=s' => \$rscriptOptions,
	   'f=s' => \$parameterFile,
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

my $command = "Rscript stoopidChromsearch.test.R $rscriptOptions";
my $autosCommand = "Rscript stoopidAutosomesearch.R $rscriptOptions";

while (1) {
    my @params = ();
    push(@params, 
	 random_uniform(1, $bounds{p1low}, $bounds{p1high}));

# loguniform
    push(@params, 
	 10**random_uniform(1, log($bounds{p2low})/$log10, log($bounds{p2high})/$log10));
    push(@params, 
	 10**random_uniform(1, log($bounds{p3low})/$log10, log($bounds{p3high})/$log10));
    push(@params, 
	 10**random_uniform(1, log($bounds{p4low})/$log10, log($bounds{p4high})/$log10));
    push(@params, 
	 10**random_uniform(1, log($bounds{p5low})/$log10, log($bounds{p5high})/$log10));


    my @results = ();

    if ($autosomes) {
	my $ans = `$autosCommand @params`;
	$ans =~ s/\s+/ /g;
	print "@params\n";
	print 'Autos '  , $ans , "\n";
    } else {
	foreach (@chroms) {
#	    my $ans = `$command chr$_ 2kb/$_.wga @params`;
#	    die "$command -c chr$_";
	    my $ans = `$command -c chr$_ @params`;
	    $ans =~ s/\s+/ /g;
	    if (length($ans) > 4) {
		print "@params\n";
		print $_ , ' '  , $ans , "\n";
	    }
	}
    }
}



