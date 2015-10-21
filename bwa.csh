/users/sturkars/samtools-1.1/
/users/sturkars/HiC-Pro/
/users/sturkars/tools/python-2.7.3/
/users/sturkars/bowtie2-2.2.3/


#DVH
## Index genome with bwa
/users/sturkars/bwa-0.7.12/bwa index dvh-genome.fna

## Align with bwa mem
/users/sturkars/bwa-0.7.12/bwa mem -R '@RG\tID:Baliga\tSM:UA3\tLB:S3' reference/dvh-genome.fna data/UA3baliga_S3HiC_R1.fastq data/UA3baliga_S3HiC_R2.fastq > dvh/UA3-dvh.sam

## Cleanup read-pair information
samtools fixmate -O bam dvh/UA3-dvh.sam dvh/UA3-dvh_fixmate.bam

## Sort them in coordinate order
samtools sort -O bam -o dvh/UA3-dvh_sorted.bam -T tmp/UA3-dvh_temp dvh/UA3-dvh_fixmate.bam

## Index bam file
samtools index dvh/UA3-dvh_sorted.bam

## Simple Stats
samtools flagstat dvh/UA3-dvh_fixmate.bam > dvh/UA3-dvh.stats

## Create fasta index file for GATK
samtools faidx reference/dvh-genome.fna

## Create dictionary file
/usr/bin/java -Xmx2g -jar /users/sturkars/picard-tools-1.139/picard.jar CreateSequenceDictionary R= reference/dvh-genome.fna O= reference/dvh-genome.fna.dict

## Realign raw gapped alignment to reduce number of INDEL miscalls
# Target interval creation
/usr/bin/java -Xmx2g -jar /users/sturkars/gatk/GenomeAnalysisTK.jar -T RealignerTargetCreator -R reference/dvh-genome.fna -I dvh/UA3-dvh_sorted.bam -o dvh/UA3-dvh.intervals

# Realignment
/usr/bin/java -Xmx4g -jar /users/sturkars/gatk/GenomeAnalysisTK.jar -T IndelRealigner -R reference/dvh-genome.fna -I dvh/UA3-dvh_sorted.bam -targetIntervals dvh/UA3-dvh.intervals -o dvh/UA3-dvh_realigned.bam

## Add readgroup information for platform
/usr/bin/java -Xmx4g -jar /users/sturkars/picard-tools-1.139/picard.jar AddOrReplaceReadGroups I=dvh/UA3-dvh_realigned.bam O=dvh/UA3-dvh_RGrealigned.bam LB="S3" PL="illumina" SM="UA3" PU="3"

# Index the bam file
samtools index dvh/UA3-dvh_RGrealigned.bam

## Detect covariates
/usr/bin/java -Xmx4g -jar /users/sturkars/gatk/GenomeAnalysisTK.jar -T BaseRecalibrator -R reference/dvh-genome.fna -knownSites dvh/UA3_varscan_variants.txt -I dvh/UA3-dvh_RGrealigned.bam -o dvh/UA3-dvh_recal.table

## Adjust quality scores
/usr/bin/java -Xmx2g -jar /users/sturkars/gatk/GenomeAnalysisTK.jar -T PrintReads -R reference/dvh-genome.fna -I dvh/UA3-dvh_RGrealigned.bam --BSQR dvh/UA3-dvh_recal.table -o dvh/UA3-dvh_recal.bam

## Mark duplicates
/usr/bin/java -Xmx2g -jar /users/sturkars/picard-tools-1.139/picard.jar MarkDuplicates.jar VALIDATION_STRINGENCY=LENIENT INPUT=dvh/UA3-dvh_recal.bam OUTPUT=dvh/UA3-dvh_library.bam

## Index the bam file
samtools index dvh/UA3-dvh_library.bam

## Produce BCF file with all locations in the genome
samtools mpileup -go dvh/UA3-dvh.bcf -f reference/dvh-genome.fna dvh/UA3-dvh_library.bam

## Reduce list of sites
bcftools call -vmO z -o dvh/UA3-dvh.vcf.gz dvh/UA3-dvh.bcf

## Prepare vcf file for querying
tabix -p vcf dvh/UA3-dvh.vcf.gz

## graph and stats
bcftools stats -F reference/dvh-genome.fna -s - dvh/UA3-dvh.vcf.gz > dvh/UA3-dvh.vcf.gz.stats
mkdir plots
plot-vcfstats -p plots/ dvh/UA3-dvh.vcf.gz.stats

## Filtering
bcftools filter -O z -o dvh/UA3-dvh.filtered-vcf.gz -s LOWQUAL -i'%QUAL>10' dvh/UA3-dvh.vcf.gz


# Mmp
## Index genome with bwa
/users/sturkars/bwa-0.7.12/bwa index mmp-genome.fna

## Align with bwa mem
/users/sturkars/bwa-0.7.12/bwa mem -R '@RG\tID:Baliga\tSM:UA3\tLB:S3' reference/mmp-genome.fna data/UA3baliga_S3HiC_R1.fastq data/UA3baliga_S3HiC_R2.fastq > UA3-mmp.sam

## Cleanup read-pair information
samtools fixmate -O bam UA3-mmp.sam UA3-mmp_fixmate.bam

## Simple Stats
samtools flagstat UA3-mmp_fixmate.bam > UA3-mmp.stats