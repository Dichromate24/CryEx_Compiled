#!/bin/bash

# This is a compiled pipeline that runs CryEx on all samples and all chromosomes automatically 

base_dir="/mnt/gtklab01/darren/CryEx_GSE227047/samples"

for SAMPLE in GSM7090418 #GSM7090417 GSM7090419 GSM7090420 GSM7090421 GSM7090422
do
    echo DOING SAMPLE $SAMPLE
    echo "--------------------- NEW ATTEMPT -----------------------" >> /mnt/gtklab01/darren/CryEx_GSE227047/samples/progress_report/${SAMPLE}_progress.txt

    for CHR in "chr1" #"chr2" "chr3" "chr4" "chr5" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chrX" "chrY" 
    do
        echo $CHR
        curr_dir=${base_dir}/${SAMPLE}/${CHR}/sub
        cd ${curr_dir}
        echo -n "STARTING:  $SAMPLE $CHR  " >> /mnt/gtklab01/darren/CryEx_GSE227047/samples/progress_report/${SAMPLE}_progress.txt
        date >> /mnt/gtklab01/darren/CryEx_GSE227047/samples/progress_report/${SAMPLE}_progress.txt

        #############################
        ### CryEx_Boundaries.sh #####
        #############################

        # This part generates an individual subfolder for each expressed exon from each chromosome of the samples
        # Each subfolder will contain the the exon's coordinate, possible 5SS and 3SS splice site coordinates and the number of reads aligned
        
        cat /mnt/gtklab01/darren/CryEx_GSE227047/all_chr_exons/${CHR}_exons.coord | while read LINE
        do
            echo $LINE
            FILENAME=$LINE
            echo $FILENAME

            # Takes in the exon coordinates and isolate the start and end coordinates
            FILE=`echo $FILENAME|awk '{print $1}'`
            echo $FILE
            DIRNAME=`echo $FILENAME|awk '{print $1"-"$2}'|sed -r 's/:/-/g'`
            echo $DIRNAME
            DOWN=`echo $FILENAME|awk '{print $1}'|sed 's/.*-//'`
            echo $DOWN
            UP=`echo $FILENAME|awk '{print $1}'|sed 's/.*://'|cut -f1 -d'-'`
            echo $UP

            # location of bam file
            BAMFILE=/mnt/gtklab01/darren/CryEx_GSE227047/all_sample_STAR/${SAMPLE}_Aligned.sortedByCoord.out.bam
            
            # Creates the Exon subdirectory using the CryEx_Coor_Exact_new.sh template
            # AXX, OXX, UXX, DXX, BAMFILE are place holders which will be replaced by the coordinates, directory information, Bam file
            sed "s/AXX/${FILE}/g" /mnt/gtklab01/darren/CryEx_GSE227047/samples/misc_scripts/CryEx_Coor_Exact_new.sh > CryEx.${FILE}.sh
            sed -i "s/OXX/${DOWN}/g" CryEx.${FILE}.sh
            sed -i "s/UXX/${UP}/g" CryEx.${FILE}.sh
            sed -i "s@DXX@${DIRNAME}@g" CryEx.${FILE}.sh
            sed -i "s@BAMFILE@${BAMFILE}@g" CryEx.${FILE}.sh

            # Runs the template script (under CryEx_Coor_Exact_new.sh) for each exon, to generate the subdirectory
            bash CryEx.${FILE}.sh

        done
        wait

        #############################
        ####### CryEx_len.sh ########
        #############################

        # this part compiles the number of 5SS and 3SS sites identified for each expressed exon in a sample
        
        touch ../output_files/CryEx_coord_len.bed
        
        for d in */
        do
            echo $d
            line=$d
            # reads in the splicing/ junction information for both splice sites for every expressed exon
            SS3_file=`printf "../sub/${line}final_junc3SS_HTE_count.txt"`
            SS5_file=`printf "../sub/${line}final_junc5SS_HTE_count.txt"`
            chr=`echo $line|cut -d'-' -f 1`
            
            # count how many 3SS and 5SS junction sites
            SS3=`echo $(cat $SS3_file|wc -l)`
            echo $SS3
            SS5=`echo $(cat $SS5_file|wc -l)`
            echo $SS5

            # the output file will contain the number of 3SS and 5SS sites each exon has
            printf "$line\t$chr\t$SS3\t$SS5\n" >> ../output_files/CryEx_coord_len.bed
        done 
        
        ######################################################
        ## Creating SJ_STAR.sh from template and running it ##
        ######################################################
        cd ../scripts
        echo NOW DOING SJ JUNC
        
        sed "s/CHROME/${CHR}/g" /mnt/gtklab01/darren/CryEx_GSE227047/samples/misc_scripts/SJ_STAR_new.sh | sed "s/GSM_SAM/${SAMPLE}/g" > SJ_STAR_new_${SAMPLE}_${CHR}.sh
        bash SJ_STAR_new_${SAMPLE}_${CHR}.sh


        ##########################################################
        ## Creating SIMPLE_JUNC.sh from template and running it ##
        ##########################################################
        
        cd ../scripts
        echo NOW DOING SIMPLE JUNC

        sed "s/CHROME/${CHR}/g" /mnt/gtklab01/darren/CryEx_GSE227047/samples/misc_scripts/SIMPLE_JUNC_new.sh | sed "s/GSM_SAM/${SAMPLE}/g" > SIMPLE_JUNC_new_${SAMPLE}_${CHR}.sh
        bash SIMPLE_JUNC_new_${SAMPLE}_${CHR}.sh

        ###################################################
        ### Creating IR.sh from template and running it ###
        ###################################################

        cd ../scripts
        # NOW DOINF IR

        sed "s/CHROME/${CHR}/g" /mnt/gtklab01/darren/CryEx_GSE227047/samples/misc_scripts/IR_new.sh | sed "s/GSM_SAM/${SAMPLE}/g" > IR_new_${SAMPLE}_${CHR}.sh
        bash IR_new_${SAMPLE}_${CHR}.sh




        echo -n "ENDED:     $SAMPLE $CHR  ">> /mnt/gtklab01/darren/CryEx_GSE227047/samples/progress_report/${SAMPLE}_progress.txt
        date >> /mnt/gtklab01/darren/CryEx_GSE227047/samples/progress_report/${SAMPLE}_progress.txt


    done


done
