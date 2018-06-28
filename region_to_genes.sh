#!/usr/bin/env bash

# use: ./region_to_genes.sh xxx.avinput
annotate_variation.pl --regionanno $1 ~/Public/ANNOVAR/annovar/humandb/ --dbtype bed --bedfile hg19_refGene_from_UCSC.sorted.bed --colsWanted 4 --buildver hg19 --outfile $(basename $1 avinput)avoutput
