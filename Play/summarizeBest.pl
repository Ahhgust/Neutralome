#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;


my $basename = 'AutosFromChrom22';
GetOptions('b=s' => \$basename);

my $i=0;
my %results;
while (<>) {
    chomp;
    my $command = $_;
    $i++;


    my $filename = "$basename.$i";
    my $data = &getDataFromFile($filename);
    if ($data eq '') {
	$filename = "$basename.redo.$i";
	$data = &getDataFromFile($filename);
	if ($data eq '') {
	    $command =~ s/>.*$//;
	    print $command, " > $basename.redo.$i\n";
#	    die;
	    next;
	}
    }

    $data =~ s/\[\d+\]//g;
    $data =~ s/^\s+//;
    my @s = split /\s+/, $data;
    foreach (@s) {
	if (substr($_, 0, 1) =~ m/[+-]/) {
	    $_ = substr($_, 1);
	}

    }
    push(@s, $i);
    $data = join ("\t", @s);
    $results{$s[0]}{$data}++;
}

my @sorted = sort {$a <=> $b} keys %results;

foreach (@sorted) {
    my @keys = keys %{$results{$_}};
    foreach my $d (@keys) {
	print STDERR $d , "\t" , $results{$_}{$d} , "\n";
    }
}

#if ($i == @results) {
#    print STDERR join("\n", @results) , "\n";
#}

sub getDataFromFile($) {
    my $filename = shift;

    if (! -f $filename) {
	return "";
    }
    if (! open IN, $filename) {
	return "";
    }

    my $ret='';
    while (<IN>) {
	chomp;
	$ret .= $_;
    }
    return $ret;
}

