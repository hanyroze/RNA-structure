#!/bin/bash

#HT 12/02/19
#script to process genomic non paired fastq datasets 
#1 make directories
#2 Run scripts

usage(){
cat << EOF
usage: $0 options

script to trim and align riboseq fastq datasets 

OPTIONS:
	-f	path to directory containing input fastq files (gzipped)
	-o	path to output dir
	-h	this help message

example usage: ./Adams_vienna_explain.sh -f . -o /export/valenfs/projects/uORFome/AdamVienna

EOF
}

while getopts ":f:o:h" opt; do
    case $opt in 
    f)
        in_dir=$OPTARG
        echo "-f input folder $OPTARG"
	;;
    o)
        out_dir=$OPTARG
        echo "-o output folder $OPTARG"
        ;;
    h)
        usage
        exit
        ;;
    ?) 
        echo "Invalid option: -$OPTARG"
        usage
        exit 1
        ;;
    esac
done

# 1. mkdir

if [ ! -d $out_dir ]; then
    mkdir $out_dir
fi
if [ ! -d ${out_dir}/windows ]; then
    mkdir ${out_dir}/windows
fi
if [ ! -d ${out_dir}/viennas ]; then
    mkdir ${out_dir}/viennas
fi
if [ ! -d ${out_dir}/csvs ]; then
    mkdir ${out_dir}/csvs
fi

# Input: transcripts in fasta (for us, only Transcripts)

# Splitter Script 1: go through all transcripts, output: new split fasta file 39 window
# Maybe try 50 window size too ? see the difference
nice -n 10 perl ./fasta_windows.pl ${in_dir} ${out_dir}/windows 39

# Vienna Wrapper Script 2: run rnaFold on all split fasta files
for file in ${out_dir}/windows/*.fa; do
   name=$(basename $file);
   prefix=${name%.fa};
   RNAfold --noPS < $file > ${out_dir}/viennas/${prefix}_vf.txt;
done

# minimum free energy collector 3: Get csv matrix of minimum free
for file in ${out_dir}/viennas/*_vf.txt; do
   nice -n 10 perl ./parse_vf.pl $file >> ${out_dir}/csvs/am.csv;
done


#Here is what I add:
#Script 4: Find best window per leader, choose it
# Think about threshold
#Script 5: Take those positions, and make a metaplot around middle point of window
# Try with 100, 100, exclude if < 100 for both directions from middle point of window.


