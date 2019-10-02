#!/usr/bin/perl -w
use strict;
#03/10/16

my $inFA=$ARGV[0];
my $outDir=$ARGV[1];
my $windowSize=$ARGV[2];

open (IN, $inFA) || die;
my $name;

while (<IN>){
	my $match=0;
	chomp;
	if (/^>(.*)$/){ #open a file for each sequence
		$name=$1;
        my $outFile=$outDir."/".$name.".fa";
        open (OUT, ">$outFile") || die;

	}else{
        my @seq = split("",$_);
        my $arrayLength=$#seq; #lats index of array (one less than size)
        my $start=0;
       
        while (($start + $windowSize-1) <= $arrayLength){ 		
            my @outSeq=();
            for ($start .. ($start+$windowSize-1)){
                push (@outSeq, $seq[$_]);
            }
      		my $outString=join("",@outSeq);
      		print OUT "$name,$outString\n";
            $start++;
        }
	} 
}
close (IN);

exit;
