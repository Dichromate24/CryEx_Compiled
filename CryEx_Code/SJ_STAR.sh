#!/bin/bash

# This part of the script will generate the splice junction information from the bam file and reads

cd ../output_files

# Obtain all spliced reads from the bam file (from the correct chromosome)
samtools view ../../../../all_sample_STAR/GSM_SAM_Aligned.sortedByCoord.out.bam | awk '($3=="CHROME")' |awk '($6~/N/)' > spliced_reads.txt

# get starting position of read, and length of initial match to calculate 5SS (similar to CryEx_Coor_Exact)
cat spliced_reads.txt | awk '{print $4}' > position.txt
cat spliced_reads.txt | awk '{print $6}' | cut -f1 -d"M" > Match.txt
paste -d '\t' position.txt Match.txt > junc5SS.txt
# The skip.txt will contain the length of skipped region, which is used to calculate the 3SS
cat spliced_reads.txt | awk '{print $6}' | cut -f1 -d"N" | sed 's/.*M//' > Skip.txt
paste -d '\t' junc5SS.txt Skip.txt > junc3SS.txt

# calculating the 5SS site and 3SS site
cat junc3SS.txt | awk '{print $1+$2-1}' > junc5SS 
cat junc3SS.txt | awk '{print $1+$2+$3-1}' > junc3SS
# Paste the 5SS and 3SS together, this will give us our splice junction
paste -d '\t' junc5SS junc3SS > junc.txt

# Now for each splice junction we count the number of reads that are supporting this junction
cat junc.txt | awk '{count[$0]++} END {for (word in count) print word, count[word]}' | sort -k1,1n > junc_count.txt
cat junc_count.txt | awk '{print $1":"$2}' > SS
paste -d '\t' SS junc_count.txt > SS_count.txt
# Fromat the file into a BED file format
cat SS_count.txt | awk '{print "chr1""\t"($2+1)"\t"$3"\t"$4}' > SS_count.tab





# This part of the script will compare the previous output with STAR aligner detected Splice junctions

# Using the STAR generated SJ.out.tab which contains the splice junction information,
# extract out the splice junction and the number of reads which are supporting the junction
cat ../../../../all_sample_STAR/GSM_SAM_SJ.out.tab |  awk '($1=="CHROME")' | awk '{print $1"\t"$2"\t"$3"\t"$7}' > SJ.tab

# with SS_count.tab file from the reads, and the SJ.tab file from STAR, we use bedtools to intersect them
# -a SJ.tab           Input SJ.tab as the input file (derived from STAR)
# -b SS_count.tab     Input SS_count.tab as the second input (derived from BAM file processing)
# -F 1      requires 100% of feature in SJ.tab overlap with a feature in SS_count.tab
# -f 1      requires 100% of feature in SS_count.tab overlap with a feature in SJ.tab
# -v        only report those entires in SJ.tab that do not have overlapping entires in SS_count.tab, it filters out junctions in SJ.tab that are also present in SS_count.tab

# We are identifying unqiue splice junctions detected by STAR but not captured by initial BAM processing and adding it in
# i presume novel or splice junctions unannotated in original BAM will be added through this step
bedtools intersect -a SJ.tab -b SS_count.tab -F 1 -f 1 -v > tmp.SJ

cat tmp.SJ SS_count.tab | awk '{print "chr1""\t"($2-1)"\t"$3"\t"$4"\t"$2":"$3}' > add_SS_count.txt
