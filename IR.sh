#!/bin/bash
cd ../output_files

# This script calculates the PSI of an intron from the splice junctions
# splice junctions does not include does generated by STAR

while read line;
do

    # extracts start position of an intron (5'SS donor of one exon)
    left_count=`echo $(echo $line|awk '{print $2}')`
    # extract end position of an intron (3'SS acceptor of another exon)
    right_count=`echo $(echo $line|awk '{print $3}')`

    # extract reads from bam file aligned with the region and save to INTRON file(as a bam file)
    samtools view -b ../../../../all_sample_STAR/GSM_SAM_Aligned.sortedByCoord.out.bam "CHROME:$left_count-$right_count" > INTRON

    # count number of reads within INTRON file
    intron_count=`echo $(samtools view INTRON|awk '(($4*1)>='"$left_count"' && ($4*1) <='"$right_count"')'|wc -l)`

    # this is the number of reads align to the exon which is associated with the current intron
    ex_count=`echo $(echo $line|awk '{print $4}')`

    # normalize intron count, adjusting for length of intron
    intron_norm=`echo "scale=2; $intron_count/($right_count - $left_count+100)"|bc`

    # calculate length of intron by subtractin starting position from ending pos and +1
    len_intron=`echo $(echo $line|awk '{print $3-$2+1}')`

    # normalise the exon count by dividing exon_count by 99
    ex_norm=`echo "scale=2; $ex_count/99"|bc`

    # calculate percent spliced in PSI (which is the normalised intron count to the sum of the normalised intron and exon counts)
    PSI=`echo "scale=2; $intron_norm/($intron_norm + $ex_norm)"|bc`

    printf "$line\t$intron_count\t$len_intron\t$PSI\n" >> ../PSI_results/CHROME_intron.PSI

done < SS_count.txt
