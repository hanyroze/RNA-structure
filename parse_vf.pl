#!/usr/bin/perl -w
use strict;
#03/10/16

my $inVF=$ARGV[0];

open (IN, $inVF) || die;

my $name;
if ($inVF =~/.*\/([^\/]+)_vf.txt/){
    $name=$1;
}elsif ($inVF=~/(.*)_vf.txt/){
    $name=$1;
}
print "$name,";

while (<IN>){
	chomp;
    
    if (/^[\.\(]/){

        if (/(\-?\d+\.\d+)\)$/){
            print "$1,";
        }        
    }
}
close (IN);

print "\n";

exit;
