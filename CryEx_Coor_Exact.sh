#!/bin/bash

### This script generates individual files that has the 5'SS and 3'SS junction reads, and junction chromosome coordinates for each expressed exon
### This script under nested with CryEx_Boundaries.sh

coord=AXX
up_coord=UXX
down_coord=OXX
datadir=DXX
bamfile=BAMFILE  

mkdir -p $datadir
cd $datadir
echo $coord

# Use samtools to obtain all the reads (from the bamfile) in the exon location 
samtools view -b ${bamfile} ${coord} > HTE

# Here we focus on spliced reads within the exon region
# Reads that were spliced will have an "N" unmatched region 
samtools view HTE | awk '($6~/N/)' > HTE_spliced_reads.txt

# In order to compute the coordinates of the 5SS junction
# we extract the starting position of the read, and concatenate the length of matched sequence length
cat HTE_spliced_reads.txt | awk '{print $4}' > position.txt
cat HTE_spliced_reads.txt | awk '{print $6}' | cut -f1 -d"M" > Match.txt
paste -d '\t' position.txt Match.txt > junc5SS_HTE.txt
# Here we compute the exact 5SS position, and count how many reads are aligned to it
cat junc5SS_HTE.txt | awk '{print $1+$2-1}' | awk '{count[$1]++} END {for (word in count) print word, count[word]}' | sort -k2,2n > junc5SS_HTE_count.txt
# then filter for only 5SS which are within 2 bases of our exon boundary (get rid of reads which were aligned to a upstream exon and spliced into current exon)
cat junc5SS_HTE_count.txt | awk '($1 >= ('${up_coord}'-2) && $1 <= ('${down_coord}'+2))' > final_junc5SS_HTE_count.txt

# Compute the exon 3SS junction coordinates in a similar way, and the number of reads aligned to it
cat HTE_spliced_reads.txt | awk '{print $6}' | cut -f1 -d"N" | sed 's/.*M//' > Skip.txt
paste -d '\t' junc5SS_HTE.txt Skip.txt > junc3SS_HTE.txt
# instead of just adding the initial read position with the length of initial match, we also add the length of skipped region
# essentially computing the position where the read spliced back in
cat junc3SS_HTE.txt | awk '{print $1+$2+$3-1}' | awk '{count[$1]++} END {for (word in count) print word, count[word]}' | sort -k2,2n > junc3SS_HTE_count.txt
# similarly we filter for 3SS within 2 bases of the exon boundary (get rid of reads which are initially aligned to the current exon and splicing out) 
cat junc3SS_HTE_count.txt | awk '($1 >= ('${up_coord}'-2) && $1 <= ('${down_coord}'+2))' > final_junc3SS_HTE_count.txt

# print out the different 5SS and 3SS location and number of reads supporting it
echo 'junc5SS'
tail final_junc5SS_HTE_count.txt
echo 'junc3SS'
tail final_junc3SS_HTE_count.txt
