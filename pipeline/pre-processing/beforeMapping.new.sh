#!/bin/bash

PIPELINE=/groups/ameres/Pooja/Projects/refiningPipeline/pipeline/





### reading in the variables passed in the command line... 

while getopts 'a: i: o: g: t:' OPTION; do
	  case "$OPTION" in
		      a)
			      avalue="$OPTARG"
			  echo "the adapter trimmed has the sequence $OPTARG"
				          ;;
	
					      i)
						      ivalue="$OPTARG"
						        echo "the input directory containing QuantSeq data is $OPTARG"
							          ;;
										
							o) 
								ovalue="$OPTARG"
								echo "the output dorectory is $OPTARG"
									;;

							
							 g) 
								 genome="$OPTARG"
								 echo "the genome is specified in $OPTARG"
								 ;;

							t) 
								threshold="$OPTARG"
								echo "the threshold to consider priming sites is $OPTARG"
								;;

					  			  ?)
										 echo "script usage: $(basename $0) [-a adapter] [-i input directory] [-o output directory] [-g genome file] [-t threshold for priming sites]" >&2
														        exit 1
															      ;;
															        esac
															done
															shift "$(($OPTIND -1))"




if [ "x" == "x$avalue" ]; then
	  echo "-a [adapter] is required"
	    exit
 fi



 if [ "x" == "x$ivalue" ]; then
	echo "-i [input directory] is required"
	exit
fi



 if [ "x" == "x$ovalue" ]; then
	         echo "-o [output directory] is required"
		         exit
fi
		     









#### i value is the input directory ... this should consist of two directories... 

	### quantseq - containing all the zipped files of the quantSeq data
	### rnaseq - containing all the rnaseq data for the condition, mapeed and indexed. 

QUANT_ALIGN=$ovalue/polyAmapping_allTimepoints
QUANT_MAP=$ovalue/polyAmapping_allTimepoints/n_100_global_a0/
mkdir -p "$QUANT_ALIGN"
mkdir -p "$QUANT_MAP"

ls "$ivalue"/quantseq/*.{gz,fastq,fq} | perl -pe "s#$ivalue/quantseq/##" > $QUANT_ALIGN/sampleInfo.txt


INPUT="$ivalue"/quantseq/

PARAMETER="$QUANT_ALIGN/sampleInfo.txt"

#set dirs

OUTDIR=$QUANT_ALIGN
LOG=$QUANT_ALIGN/logs/
mkdir -p $LOG




module load cutadapt/1.9.1-foss-2017a-python-2.7.13

 #arrayfile=`ls $INPUT | awk -v line=$SLURM_ARRAY_TASK_ID '{if (NR == line) print $0}'`
#cutadapt -a $avalue -o "$OUTDIR"/"$arrayfile"_trimmed.fastq  --trim-n $INPUT/"$arrayfile" 2>"$LOG"/stderr_"$index".txt 1>>"$LOG"/stdo_"$index".txt





#init files
touch "$OUTDIR"/logFile.txt ## initialise a log file
touch "$OUTDIR"/readStats.txt ## initialise a file to store your read stats at different processing steps
touch "$OUTDIR"/processingSteps.txt ### just another file to list the processing steps 

touch "$LOG"/stderr_"$index".txt
touch "$LOG"/stdo_"$index".txt


date >"$OUTDIR"/logFile.txt
 date >>"$LOG"/stderr_"$index".txt
 date >>"$LOG"/stdo_"$index".txt

#################
#ADAPTER trimming
#################

module load python/2.7.13-foss-2017a
module load cutadapt/1.9.1-foss-2017a-python-2.7.13
module load fastx-toolkit/0.0.14-foss-2017a
module load fastqc/0.11.5-java-1.8.0_121
module load anaconda2/5.1.0
module load bedtools/2.27.1-foss-2017a
module load r/3.4.1-foss-2017a-x11-20170314


cutadapt_version=$(cutadapt --version)
echo running cutadapt version "$cutadapt_version" >> "$LOG"/stdo_"$index".txt

## the step below trims the adapter : you should replace this based on your library prep method

while read index; do
	 echo "processing $index"
  
	  ### trimming adapter and adding the relavant information. 

  	#######cutadapt -a $avalue -o "$OUTDIR"/"$index"_trimmed.fastq  --trim-n $INPUT/"$index" 2>"$LOG"/stderr_"$index".txt 1>>"$LOG"/stdo_"$index".txt
  	
	adapterTrimming=$(cat  "$OUTDIR"/"$index"_trimmed.fastq | echo $((`wc -l`/4)))
  	echo the number of reads after adapter trimming is "$adapterTrimming" >>"$LOG"/stdo_"$index".txt
	
	#### just adding an N to empty lines.... 
	
	sed 's/^$/N/' "$OUTDIR"/"$index"_trimmed.fastq > "$OUTDIR"/"$index"_trimmed_emptyRemoved.fastq



	### QC  
	
	#module load fastqc/0.11.5-java-1.8.0_121
	#fastqc "$OUTDIR"/"$index"_trimmed.fastq 
	
	### trimming from the 5' end....

	echo fastx_trimmer to remove 12 nts from the 5 end of the reads >>"$LOG"/stdo_"$index".txt  
	echo version of fastx_trimmer used is Version: `which fastx_trimmer` >>"$LOG"/stdo_"$index".txt
	###### fastx_trimmer -Q33 -f 13 -m 1 -i "$OUTDIR"/"$index"_trimmed_emptyRemoved.fastq   > "$OUTDIR"/"$index"_5primetrimmed_trimmed.fastq 2>"$LOG"/stderr_"$index".txt 
	fivePrimeTrimmed=$(cat  "$OUTDIR"/"$index"_5primetrimmed_trimmed.fastq | echo $((`wc -l`/4)))                                                                                    
	echo number of reads after 5 trimmig "$fivePrimeTrimmed" >>"$LOG"/stdo_"$index".txt
	echo "retaining reads that have >=5 As the the 3 end" >>"$LOG"/stdo_"$index".txt


	#### removing A's from the 3' end

	#AAAAA$ could theoretically also happen in PHRED score
	egrep -A2 -B1 'AAAAA$' "$OUTDIR"/"$index"_5primetrimmed_trimmed.fastq  | sed '/^--$/d' > "$OUTDIR"/"$index"_5primetrimmed_trimmed_sizefiltered_polyAreads.fastq  2>"$LOG"/stderr_"$index".txt
	polyAreads=$(cat "$OUTDIR"/"$index"_5primetrimmed_trimmed_sizefiltered_polyAreads.fastq | echo $((`wc -l`/4)))
	echo readd withs polyA "$polyAreads"  >>"$LOG"/stdo_"$index".txt
	echo remove the polyAs at the end of polyA reads >>"$LOG"/stdo_"$index".txt
	
	#echo this is done using cutadapat. At this step we also perform a size filteirng of minimum 18 nucleotides >>"$LOG"/stdo_"$index".txt
	#cut super long polyA to avoid internal polyA cutting of 5 or more As
	###### cutadapt --no-indels -m 18 -e 0 -a "A{1000}"  -o "$OUTDIR"/"$index"_5primetrimmed_trimmed_sizefiltered_polyAreads_polyAremoved.fastq  "$OUTDIR"/"$index"_5primetrimmed_trimmed_sizefiltered_polyAreads.fastq 1>>"$LOG"/stdo_"$index".txt 2>"$LOG"/stderr_"$index".txt
	finalReads=$(cat "$OUTDIR"/"$index"_5primetrimmed_trimmed_sizefiltered_polyAreads_polyAremoved.fastq | echo $((`wc -l`/4)))
	echo final reads after size filtering "$finalReads"  >>"$LOG"/stdo_"$index".txt

	################
	# write stats
	#################

	touch "$LOG"/preProcessingNumbers_"$index".txt

	echo initialFile:"$initialrReads" >>"$LOG"/preProcessingNumbers_"$index".txt
	echo adapterTrimmed:"$adapterTrimming" >>"$LOG"/preProcessingNumbers_"$index".txt
	echo fivePrimeTrimming:"$fivePrimeTrimmed" >>"$LOG"/preProcessingNumbers_"$index".txt
	echo polyAcontaining:"$polyAreads" >>"$LOG"/preProcessingNumbers_"$index".txt
	echo finalFile:"$finalReads" >>"$LOG"/preProcessingNumbers_"$index".txt
	echo the pre processing before mapping has been completed.  >>"$LOG"/stdo_"$index".txt


	echo initialFile:"$initialrReads"  >>"$LOG"/stdo_"$index".txt
	echo adapterTrimmed:"$adapterTrimming"  >>"$LOG"/stdo_"$index".txt
	echo fivePrimeTrimming:"$fivePrimeTrimmed"  >>"$LOG"/stdo_"$index".txt
	echo polyAcontaining:"$polyAreads"  >>"$LOG"/stdo_"$index".txt
	echo finalFile:"$finalReads"  >>"$LOG"/stdo_"$index".txt


	#################
	# getting the read length distribution 
	#################

	awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}' "$OUTDIR"/"$index"_5primetrimmed_trimmed_sizefiltered_polyAreads_polyAremoved.fastq >"$OUTDIR"/"$index"_lengthDistribution.txt
	
	echo this is done
	

	###### mapping the data 	
		
		export PYTHONNOUSERSITE=1
		 source activate slamdunk0.3.3

		#slamdunk map -r $genome -o $QUANT_MAP/ -n 100 -5 0 -a 0 -e "$QUANT_ALIGN"/"$index"_5primetrimmed_trimmed_sizefiltered_polyAreads_polyAremoved.fastq
		#slamdunk filter -o $QUANT_MAP -mq 0 -mi 0.95  $QUANT_MAP/"$index"_5primetrimmed_trimmed_sizefiltered_polyAreads_polyAremoved_slamdunk_mapped.bam         

		source deactivate
	
		
		 #bedtools bamtobed -i $QUANT_MAP/"$index"_5primetrimmed_trimmed_sizefiltered_polyAreads_polyAremoved_slamdunk_mapped.bam   > $QUANT_MAP/"$index"_5primetrimmed_trimmed_sizefiltered_polyAreads_polyAremoved_slamdunk_mapped.bam_bamTobed.bed 

		 #awk -vFS="\t" '$6 == "-"' $QUANT_MAP/"$index"_5primetrimmed_trimmed_sizefiltered_polyAreads_polyAremoved_slamdunk_mapped.bam_bamTobed.bed > $QUANT_MAP/"$index"_5primetrimmed_trimmed_sizefiltered_polyAreads_polyAremoved_slamdunk_mapped.bam_bamTobed_minusStrand.bed
		  
		  # awk -vFS="\t" '$6 == "+"' $QUANT_MAP/"$index"_5primetrimmed_trimmed_sizefiltered_polyAreads_polyAremoved_slamdunk_mapped.bam_bamTobed.bed > $QUANT_MAP/"$index"_5primetrimmed_trimmed_sizefiltered_polyAreads_polyAremoved_slamdunk_mapped.bam_bamTobed_plusStrand.bed




  done <"$QUANT_ALIGN/sampleInfo.txt"





########################## Once this is done. I need to find the priming site ##################################



#################### identification of priming sites... 



cat $QUANT_MAP/*_bamTobed_minusStrand.bed | awk -vFS="\t" '{print $1,$2}'  | sort | uniq -c | perl -pe 's#^\s+##'  > $QUANT_MAP/polyAreads_polyAremoved_pooled_slamdunk_mapped_filtered_bamTobed_minusStrand_countsUnique.bed
cat $QUANT_MAP/*_bamTobed_plusStrand.bed | awk  -vFS="\t" '{print $1,$3}' | sort | uniq -c | perl -pe 's#^\s+##' > $QUANT_MAP/polyAreads_polyAremoved_pooled_slamdunk_mapped_filtered_bamTobed_plusStrand_countsUnique.bed


onPositive=$(wc -l  "$QUANT_MAP"/polyAreads_polyAremoved_pooled_slamdunk_mapped_filtered_bamTobed_plusStrand_countsUnique.bed)
onNegative=$(wc -l  "$QUANT_MAP"/polyAreads_polyAremoved_pooled_slamdunk_mapped_filtered_bamTobed_minusStrand_countsUnique.bed)

echo the numeber of unique positions of the genome that are covered by reads on the plus strand is "$onPositive" and on the minus strand is "$onNegative" 





awk -vT=$threshold '{ if ($1 >=T) print  }' "$QUANT_MAP"/polyAreads_polyAremoved_pooled_slamdunk_mapped_filtered_bamTobed_minusStrand_countsUnique.bed >"$QUANT_MAP"/polyAreads_polyAremoved_pooled_slamdunk_mapped_filtered_bamTobed_minusStrand_countsUnique_greaterThan"$threshold".bed
awk -vT=$threshold '{ if ($1 >=T) print  }' "$QUANT_MAP"/polyAreads_polyAremoved_pooled_slamdunk_mapped_filtered_bamTobed_plusStrand_countsUnique.bed >"$QUANT_MAP"/polyAreads_polyAremoved_pooled_slamdunk_mapped_filtered_bamTobed_plusStrand_countsUnique_greaterThan"$threshold".bed


awk -vOFS="\t" '{print $2, $3-1, $3, $1}' $QUANT_MAP"/polyAreads_polyAremoved_pooled_slamdunk_mapped_filtered_bamTobed_plusStrand_countsUnique_greaterThan"$threshold".bed" | sort -k1,1 -k2,2n  > "$QUANT_MAP"/polyAreads_polyAremoved_pooled_slamdunk_mapped_filtered_bamTobed_plusStrand_countsUnique_greaterThan"$threshold".bed_sorted.bed_changedCoordinates.bed
awk -vOFS="\t" '{print $2, $3, $3+1, $1}' $QUANT_MAP"/polyAreads_polyAremoved_pooled_slamdunk_mapped_filtered_bamTobed_minusStrand_countsUnique_greaterThan"$threshold".bed" | sort -k1,1 -k2,2n  > "$QUANT_MAP"/polyAreads_polyAremoved_pooled_slamdunk_mapped_filtered_bamTobed_minusStrand_countsUnique_greaterThan"$threshold".bed_sorted.bed_changedCoordinates.bed




bedtools merge -d 0 -i "$QUANT_MAP"/polyAreads_polyAremoved_pooled_slamdunk_mapped_filtered_bamTobed_minusStrand_countsUnique_greaterThan"$threshold".bed_sorted.bed_changedCoordinates.bed -c 4 -o count,collapse >   "$QUANT_MAP"/polyAreads_polyAremoved_pooled_slamdunk_mapped_filtered_bamTobed_minusStrand_countsUnique_greaterThan"$threshold".bed_sorted_merged.bed

bedtools merge  -d 0  -i "$QUANT_MAP"/polyAreads_polyAremoved_pooled_slamdunk_mapped_filtered_bamTobed_plusStrand_countsUnique_greaterThan"$threshold".bed_sorted.bed_changedCoordinates.bed -c 4 -o count,collapse > "$QUANT_MAP"/polyAreads_polyAremoved_pooled_slamdunk_mapped_filtered_bamTobed_plusStrand_countsUnique_greaterThan"$threshold".bed_sorted_merged.bed 
 

######################################
###3 getting the sequences +/- 60 nts around the priming sites... 
######################################


Rscript  --vanilla $PIPELINE/primingSites/overlappingPolyApeaks.R "$QUANT_MAP"/polyAreads_polyAremoved_pooled_slamdunk_mapped_filtered_bamTobed_plusStrand_countsUnique_greaterThan"$threshold".bed_sorted_merged.bed "$QUANT_MAP"/polyAreads_polyAremoved_pooled_slamdunk_mapped_filtered_bamTobed_minusStrand_countsUnique_greaterThan"$threshold".bed_sorted_merged.bed "$QUANT_MAP"/peaks_"$threshold"_120bps.bed 



bedtools getfasta -s -fi $genome -bed "$QUANT_MAP"/peaks_"$threshold"_120bps.bed > "$QUANT_MAP"/peaks_"$threshold"_120bps.fa

#### getting sequences in the 120 nucleotide window

rmd="$PIPELINE/OverlappingPrimingSitesWithAnnotations/sequencesForNucleotideProfile.R"
Rscript --vanilla -e "InPath='$QUANT_MAP';threshold="$threshold";source('$rmd')"

#####################
#### creating nucleotie profile plots...



#############3

















































