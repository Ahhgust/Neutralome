#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

my $k=10;
my $maxRecords=0; # the max number of records to consider (0 is treated as +inf)
GetOptions('k=i' => \$k,
	   'm=i' => \$maxRecords
);

my %big = ();
my $total=0;
while (<>) {
    chomp;
    my @s = split /\s+/, $_;
    if (@s > 2 && $s[1] =~ m/\[/) {

	my $chrom = shift @s;
	shift @s;
	next if ! @s || $s[0] < 1; 
	
	foreach (@s) {
	    die "@s \n" if $_ =~ m/max/i;
	    $_ = -$_ if $_ !~ m/\[/ && $_ < 0.0;
	}
	$total++;
	my $str = join(" ", @s);
	if (! exists $big{$chrom}{$str}) {
	    $big{$chrom}{$str} = [$s[0], 1];
	} else {
	    $big{$chrom}{$str}[1]++; # count the number of times we've seen this point in the parameter space
	}
	last if $total ==$maxRecords;
    }
}


my @chroms = sort keys %big;
print "Total records: $total\n";
foreach my $c (@chroms) {
    my $hash = $big{$c};
    my @keys = sort {$hash->{$a}->[0] <=> $hash->{$b}->[0]} keys %$hash; # sort by the sum of squares
    my $i=0;
    foreach my $record (@keys) {
	print "$c $record " , $hash->{$record}->[1] , "\n";
	$i++;
	last if $i >= $k;
    }
}

