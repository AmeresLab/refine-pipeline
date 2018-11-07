# Annotation pipeline

This is a set of scripts used to refine and extend existing mRNA 3' end annotations using Quant seq 3' end sequencing 	 (https://www.nature.com/articles/nmeth.f.376) and RNAseq. 

## Installation

Clone from github
git clone https://github.com/poojabhat1690/refine-pipeline.git
cd pipeline/pre-processing/

The following dependencies that have to be installed. 
			1. cutadapt
			2. bedtools
			3. python
			4. R

## Quickstart

Please make sure to download annotations and format them before running this script. 
The script required to run the whole pipeline is beforeMapping.new.sh

beforeMapping.new.sh -a [adapter] -i [input directory] -o [output directory] -g [genome file] -t [threshold for priming sites]
-u [ucscDir] -e [ensemblDir] -m [mode rnaseq p/s/S] -c [condition]
 
 
 -a 3' adapter sequences that has to be removed using cutadapt
 -i input directory containing two folders named - quantseq, rnaseq
                   quantseq: contains *.fastq, *.fq.gz, *fq  files. 
                   rnaseq: mapped, sorted and indexed bam files. 
 -o Output directory
 -t threshold of the number of reads to define a priming site.
 -u ucsc directory containing annotations downloaded from ucsc table browser. 
 -e ensembl directory containing ensembl annotations ontained from biomart. 
 -m mode of counting for RNAseq coverage, derived from bedtools multicov (s: counting on the same strand, 
            p: paired end reads, S: counting on the opposite strand
 -c condition of sample (example: timepoint or organism)
                
           
 ## Prerequisites: 
 
 
 
 1. Annotations from refSeq, downloaded from the UCSC table browser. 
 
 
 As an example, taking the mm10 annotation (the refSeq annotations were downloaded manually and processed) : File name: getAnnotations.Rmd
    for 3' UTR annotations from UCSC genome browser (refSeq_mm10_3primeUTR.bed) 
    
            1. Form UCSC table browser, select:
                a. clade = mammal 
                b. genome = mouse
                c. assembly = Dec. 2011 (GRCm38/mm10)
                d. group = genes and gene prediction
                e. track = refSeq genes
                f. table = refGene
              2. select bed format and 3' UTR. 
              
       for intron annotations from UCSC table browser (refSeq_mm10_intronAnnotations.bed)
     
          	1. Form UCSC table browser, select:
              a. clade = mammal 
              b. genome = mouse
              c. assembly = Dec. 2011 (GRCm38/mm10)
              d. group = genes and gene prediction
              e. track = refSeq genes
              f. table = refGene
            2. select bed format and intron. 

      for exon annotations from UCSC table browser (refSeq_mm10_exonAnnotations.bed)
     
            1. Form UCSC table browser, select:
              a. clade = mammal 
              b. genome = mouse
              c. assembly = Dec. 2011 (GRCm38/mm10)
              d. group = genes and gene prediction
              e. track = refSeq genes
              f. table = refGene
            2. select bed format and exon. 

      refFlat annotations from the UCSC genome browser 
     
	
            1. Form UCSC table browser, select:
              a. clade = mammal 
              b. genome = mouse
              c. assembly = Dec. 2011 (GRCm38/mm10)
              d. group = genes and gene prediction
              e. track = refSeq genes
              f. table = refFlat
            2. select bed format and 3' UTR. 
further processing has been done in the script : pipeline/pre-processing/getAnnotations.Rmd
     
1. assign gene names of refSeq annotations from refFlat annotations - this is done based on the chromosome, start and end positions, only retain annotations of main chromosomes (1:19,X,y)
2. separate refSeq mRNA annotations (id : "NM...") and non-coding RNA annotations (id : NR...) - refSeq_mrna_utrsPresent.txt,  refSeq_ncrna_utrsPresent.txt
3. Get the transcript annotations and check which transcripts do not have an annotated 3' UTR. - refSeqGenesWithoutUTRs.txt
	


 2. Annotations from ENSEMBL retreived from biomaRt
 
Ensembl annotations were retrieved from biomart, using RStudio using the script pipeline/pre-processing/getAnnotations.Rmd
processing ensembl annotations 
1. get only protein coding genes based on transcript biotype : allTranscripts_proteinCoding.txtonly the main chromosomes are retained.
2. get protein coding genes with annotated 3' UTRs : proteinCoding_annotatedUTRs.txt
3. protein coding genes with un annotated 3' UTRs : proteinCoding_UnannotatedUTRs.txt
	

 

## Output description

The ouput and all intermediate files are organized in the follwing folders:

polyAmapping_allTimepoints : contains raw files of the pre-processing steps, including adapter trimmed, poly(A) tail filtered fastq files, mapped files, priming sites, priming sites overalapping with different annotations (ENSEMBL 3' UTRs, refSeq 3' UTRs, ENSEMBL introns, ENSEMBL exons, nonOverlapping).


PASplots - contains nucleotide profiles for priming sites overlapping with annotations, sparated by presence or absence of the poly A signal (PAS) and separated by downstream genomic A content.

coverage : contains list of priming sites (filterd for low genomic A content), that overlap with un-annotated regions supported by RNAseq.

final90percent :

1. ends_greater90percent_intergenic_n100 : contains the high condience mRNA 3' ends

2. allAnnotations.bed : 250nt counting windows (overlapping counting windows merged), used to count quantSeq reads. 

3. countingWindows_transcriptionalOutput.bed : genomic loci to filter multimappers using SLAMdunk. These include all counting windows + all 3' UTRs + all extended counting windows. 

4. onlyIntergenic_90percent_n100: list of the intergenic counting windows created using presence of continuous RNAseq signal.


## Applications:
https://www.nature.com/articles/nmeth.4435
