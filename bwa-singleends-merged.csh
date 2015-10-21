#!/bin/bash
## We combine paired end reads into single file and then align to dvh and mmp genomes by using bwa-gatk pipeline
RAW_DATA="/Volumes/omics4tb/HiC/data/UA3baliga_S3HiC_merged.fastq"
DATA_DIR="/Volumes/omics4tb/HiC/data/UA3-bwa"

## Combine both read files into single fastq
#cat $DATA_DIR/UA3baliga_S3HiC_R1.fastq $DATA_DIR/UA3baliga_S3HiC_R2.fastq > $DATA_DIR/UA3baliga_S3HiC_merged.fastq

#DVH
## Index genome with bwa
echo -e "\033[34m Indexing genome...\033[0m"
/Users/sturkars/apps/bwa-0.7.12/bwa index /Users/sturkars/apps/hicup_reference/dvh-genome.fasta

## Align with bwa mem
echo -e "\033[34m Performing alignment with bwa...\033[0m"
/Users/sturkars/apps/bwa-0.7.12/bwa mem -M -R  '@RG\tID:Baliga\tSM:UA3\tPL:illumina\tLB:S3\tPU:unit3' -p /Users/sturkars/apps/hicup_reference/dvh-genome.fasta $RAW_DATA > $DATA_DIR/UA3_bwa_merged_aligned.sam

## Sort SAM file and convert to BAM
echo -e "\033[34m Sorting SAM file and converting to BAM...\033[0m"
/usr/bin/java -Xmx2g -jar /Users/sturkars/apps/picard-tools-1.140/picard.jar SortSam INPUT=$DATA_DIR/UA3_bwa_merged_aligned.sam OUTPUT=$DATA_DIR/UA3_bwa_merged_sorted.bam SORT_ORDER=coordinate

## Mark duplicates
echo -e "\033[34m Marking Duplicates...\033[0m"
/usr/bin/java -Xmx2g -jar /Users/sturkars/apps/picard-tools-1.140/picard.jar MarkDuplicates INPUT=$DATA_DIR/UA3_bwa_merged_sorted.bam OUTPUT=$DATA_DIR/UA3_bwa_merged_dedup.bam METRICS_FILE=$DATA_DIR/UA3_bwa_metrics.txt

## Index BAM file
echo -e "\033[34m Indexing BAM file...\033[0m"
/usr/bin/java -Xmx2g -jar /Users/sturkars/apps/picard-tools-1.140/picard.jar BuildBamIndex INPUT=$DATA_DIR/UA3_bwa_merged_dedup.bam

# Summary Metric for deduplicated BAM
echo -e "\033[34m Calculating summary metrics for dedup BAM file...\033[0m"
/usr/bin/java -Xmx2g -jar /Users/sturkars/apps/picard-tools-1.140/picard.jar CollectAlignmentSummaryMetrics R=/Users/sturkars/apps/hicup_reference/dvh-genome.fasta I=$DATA_DIR/UA3_bwa_merged_dedup.bam O=$DATA_DIR/UA3_bwa_merged_dedup_metrics.txt 

### Realign raw gapped alignment to reduce number of INDEL miscalls
# Create Genome dictionary
echo -e "\033[34m Creating Genome Dictionary...\033[0m"
/usr/bin/java -Xmx2g -jar /Users/sturkars/apps/picard-tools-1.140/picard.jar CreateSequenceDictionary REFERENCE=/Users/sturkars/apps/hicup_reference/dvh-genome.fasta OUTPUT=/Users/sturkars/apps/hicup_reference/dvh-genome.dict

# Create Genome Fasta index
echo -e "\033[34m Creating Genome Fasta index...\033[0m"
/Users/sturkars/apps/samtools-1.2/samtools faidx /Users/sturkars/apps/hicup_reference/dvh-genome.fasta

# Target interval creation
echo -e "\033[34m Creating Target Intervals...\033[0m"
/usr/bin/java -Xmx2g -jar /Users/sturkars/apps/gatk/GenomeAnalysisTK.jar -T RealignerTargetCreator -R /Users/sturkars/apps/hicup_reference/dvh-genome.fasta -I $DATA_DIR/UA3_bwa_merged_dedup.bam -o $DATA_DIR/UA3_bwa_merged_realingment_targets.list

# Realignment
echo -e "\033[34m Performing Realignment...\033[0m"
/usr/bin/java -Xmx4g -jar /Users/sturkars/apps/gatk/GenomeAnalysisTK.jar -T IndelRealigner -R /Users/sturkars/apps/hicup_reference/dvh-genome.fasta -I $DATA_DIR/UA3_bwa_merged_dedup.bam -targetIntervals $DATA_DIR/UA3_bwa_merged_realingment_targets.list -o $DATA_DIR/UA3_bwa_merged_realigned.bam

### recalibration of base quality scores
## Detect covariates 1st Pass
echo -e "\033[34m Detecting covariates 1st Pass...\033[0m"
/usr/bin/java -Xmx4g -jar /Users/sturkars/apps/gatk/GenomeAnalysisTK.jar -T BaseRecalibrator -R /Users/sturkars/apps/hicup_reference/dvh-genome.fasta -knownSites $DATA_DIR/UA3_varscan_variants.txt -I $DATA_DIR/UA3_bwa_merged_realigned.bam -o $DATA_DIR/UA3_bwa_merged_recal_data.table

## Detect covariates 2nd Pass
echo -e "\033[34m Detecting covariates 2nd Pass...\033[0m"
/usr/bin/java -Xmx4g -jar /Users/sturkars/apps/gatk/GenomeAnalysisTK.jar -T BaseRecalibrator -R /Users/sturkars/apps/hicup_reference/dvh-genome.fasta -knownSites $DATA_DIR/UA3_varscan_variants.txt -BQSR $DATA_DIR/UA3_bwa_merged_recal_data.table -I $DATA_DIR/UA3_bwa_merged_realigned.bam -o $DATA_DIR/UA3_bwa_merged_post_recal_data.table

# Generate plots
echo -e "\033[34m Generating recalibration plots...\033[0m"
/usr/bin/java -Xmx4g -jar /Users/sturkars/apps/gatk/GenomeAnalysisTK.jar -T AnalyzeCovariates -R /Users/sturkars/apps/hicup_reference/dvh-genome.fasta -before $DATA_DIR/UA3_bwa_merged_recal_data.table -after $DATA_DIR/UA3_bwa_merged_post_recal_data.table -plots $DATA_DIR/UA3_bwa_merged_recalibrationPlots.pdf

# Apply calibration to sequence data
echo -e "\033[34m Applying calibration to sequence data...\033[0m"
/usr/bin/java -Xmx2g -jar /Users/sturkars/apps/gatk/GenomeAnalysisTK.jar -T PrintReads -R /Users/sturkars/apps/hicup_reference/dvh-genome.fasta -I $DATA_DIR/UA3_bwa_merged_realigned.bam -BQSR $DATA_DIR/UA3_bwa_merged_recal_data.table -o $DATA_DIR/UA3_bwa_merged_recal.bam

# Summary Metric for final bam
echo -e "\033[34m Calculating summary metrics for recalibrated BAM file...\033[0m"
/usr/bin/java -Xmx2g -jar /Users/sturkars/apps/picard-tools-1.140/picard.jar CollectAlignmentSummaryMetrics R=/Users/sturkars/apps/hicup_reference/dvh-genome.fasta I=$DATA_DIR/UA3_bwa_merged_recal.bam O=$DATA_DIR/UA3_bwa_merged_recal_summary.txt

# WGS Metrics for final bam
echo -e "\033[34m Calculating WGS summary metrics for recalibrated BAM file...\033[0m"
/usr/bin/java -Xmx2g -jar /Users/sturkars/apps/picard-tools-1.140/picard.jar CollectWgsMetrics R=/Users/sturkars/apps/hicup_reference/dvh-genome.fasta I=$DATA_DIR/UA3_bwa_merged_recal.bam O=$DATA_DIR/UA3_bwa_merged_recal_WGSsummary.txt


# Call variants
echo -e "\033[34m Calling variants...\033[0m"
/usr/bin/java -Xmx2g -jar /Users/sturkars/apps/gatk/GenomeAnalysisTK.jar -T HaplotypeCaller -R /Users/sturkars/apps/hicup_reference/dvh-genome.fasta -I $DATA_DIR/UA3_bwa_merged_recal.bam --genotyping_mode DISCOVERY -stand_emit_conf 10 -stand_call_conf 30 -o $DATA_DIR/UA3_bwa_merged_raw_variants.vcf