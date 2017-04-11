#!/usr/bin/env perl
use strict;
use warnings;

my %hashy = ();
my $printAt=2;
while (<>) {

    if ($_ =~ m/^chr/) {
	chomp;
	my @s = split /\s+/, $_;
	my $key;
	if ($s[0] !~ m/:/) {
	    $key = "$s[0]:$s[1]-$s[2]";
	} else {
	    $key = $s[0];
	}
	push(@{$hashy{$key}}, \@s);
	if (@{$hashy{$key}} == $printAt) {
	    
	    foreach my $val (@{$hashy{$key}}) {
		
		print join(" ", @$val) , " ";
	    }

	    print "\n";
	    delete $hashy{$s[0]};
	}
    }
}

