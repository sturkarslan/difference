### this script is a pipeline for the analysis of mutations from raw reads

import os,sys,shutil

def directoryStructureInitializer():

    #! create a tmp directory and an output directory for the results files
    if os.path.exists(tmpDir):
        shutil.rmtree(tmpDir)
    os.mkdir(tmpDir)

    #if not os.path.exists(outputDir):
    #    os.mkdir(outputDir)

    #! make n directories for the reads files
    files=os.listdir(dataPath)
    firstPairs=[]
    secondPairs=[]
    for file in files:
        if file.count('_1'):

            firstPair=dataPath+'/'+file
            secondPair=dataPath+'/'+file.replace('_1','_2')

            fileLabel=file.split('_1')[0].split('/')[-1]
            targetDir=tmpDir+'/'+fileLabel

            os.mkdir(targetDir)
            shutil.copy(firstPair,targetDir)
            shutil.copy(secondPair,targetDir+'/mates.fastq')

            firstPairs.append(firstPair)
            secondPairs.append(secondPair)

            ### making final output directories
            tag=fileLabel.split('_')[-1]
            finalDir=outputDir+'/'+tag
            if not os.path.exists(finalDir):
                os.mkdir(finalDir)

    if len(firstPairs) != len(secondPairs):
        print 'error. read paired-file missing. Exiting...'

    #! copy the genome fasta file
    shutil.copy(genomePath,tmpDir)
    shutil.copy(variantsPath,tmpDir)

    return firstPairs

def genomeWorker():

    tmpGenomePath=tmpDir+'/'+genomePath.split('/')[-1]
    # indexing genomes
    cmd1=bwaPath+' index '+tmpDir+'/'+genomePath.split('/')[-1]
    cmd2=samtoolsPath+' faidx '+tmpDir+'/'+genomePath.split('/')[-1]
    cmd3='%s -Xmx2g -jar %s CreateSequenceDictionary R=%s O=%s.dict'%(javaPath,piccardPath,tmpGenomePath,tmpGenomePath.replace('.fasta',''))

    os.system(cmd1)
    os.system(cmd2)
    os.system(cmd3)
    
    return None

def pipeline(firstPairs):

    print '\trunning the pipeline...'
    print

    # 1. prepare genome indexes
    print '\tcreating genome indexes and dictionaries...'
    print
    genomeWorker()
    
    # 2. run the data cleanup (sequential)
    print
    print '\trunning sequential analysis...'
    for element in firstPairs:
        sequentialVariantCaller(element)

    # 3. run the joint genotyping (ensemble) GenotypeGVCFs
    print
    print '\t running the joint genotyping...'
    print
    gVCF_files=[]
    tmpGenomePath=tmpDir+'/'+genomePath.split('/')[-1]
    
    for element in firstPairs:
        pathID=element.split('_1')[0].split('/')[-1]
        workingFile=tmpDir+'/'+pathID+'/file'
        
        theFile='%s_output.raw.snps.indels.g.vcf'%workingFile
        gVCF_files.append(theFile)
        
    allFiles=' --variant '.join(gVCF_files)
    
    cmd='%s -Xmx2g -jar %s -T GenotypeGVCFs -R %s --variant %s -o %s/ensembleOutput.vcf'%(javaPath,gatkPath,tmpGenomePath,allFiles,tmpDir)
    os.system(cmd)

    # 4. run the joint variant filtering (ensemble) VQSR: applying hard filters, no true set available for yeast
    print
    print '\t running the joint variant filtering (hard filters because unavailable variant databases for yeast)...'
    print

    # 4.1. extract the SNPs from the call set
    cmd1='%s -Xmx2g -jar %s -T SelectVariants -R %s -V %s/ensembleOutput.vcf -selectType SNP -o %s/raw_snps.vcf'%(javaPath,gatkPath,tmpGenomePath,tmpDir,tmpDir)
    
    # 4.2. apply the filter to the SNPs
    cmd2='%s -Xmx2g -jar %s -T VariantFiltration -R %s -V %s/raw_snps.vcf --filterName "ADL_defaultFilter_SNP" --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || InbreedingCoeff < -0.8 || SOR > 4.0" -o %s/filtered_snps.vcf'%(javaPath,gatkPath,tmpGenomePath,tmpDir,tmpDir) # control for || SOR > 4.0?

    # 4.3. extract the indels
    cmd3='%s -Xmx2g -jar %s -T SelectVariants -R %s -V %s/ensembleOutput.vcf -selectType INDEL -o %s/raw_indels.vcf'%(javaPath,gatkPath,tmpGenomePath,tmpDir,tmpDir)

    # 4.4. apply the filter to the indels
    cmd4='%s -Xmx2g -jar %s -T VariantFiltration -R %s -V %s/raw_indels.vcf --filterName "ADL_defaultFilter_INDEL" --filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0 || InbreedingCoeff < -0.8 || SOR > 10.0" -o %s/filtered_indels.vcf'%(javaPath,gatkPath,tmpGenomePath,tmpDir,tmpDir) 
    
    print cmd1
    print
    os.system(cmd1)
    print
    print cmd2
    print
    os.system(cmd2)
    print
    print cmd3
    print
    os.system(cmd3)
    print
    print cmd4
    print
    os.system(cmd4)
    print

    return None

def sampleSpecificAnalysis():

    '''
    this function breaks the ensemble into sample specific files.
    '''

    tmpGenomePath=tmpDir+'/'+genomePath.split('/')[-1]
    for pair in firstPairs:
        x=pair.split('_1')[0]
        theSampleName=x.split('_')[-1]
    
        # 1. calling sample specific SNPs
        cmd1='%s -Xmx2g -jar %s -T SelectVariants -R %s -V %s/filtered_snps.vcf -o %s/%s_specific_snps.vcf -sn %s'%(javaPath,gatkPath,tmpGenomePath,tmpDir,tmpDir,theSampleName,theSampleName)

        # 2. calling sample specific indels
        cmd2='%s -Xmx2g -jar %s -T SelectVariants -R %s -V %s/filtered_indels.vcf -o %s/%s_specific_indels.vcf -sn %s'%(javaPath,gatkPath,tmpGenomePath,tmpDir,tmpDir,theSampleName,theSampleName)

        print cmd1
        print
        os.system(cmd1)
        print
        print cmd2
        os.system(cmd2)
        print

        # 3. post hoc analysis
        sampleOutputDir=outputDir+'/%s'%theSampleName

        # 3.1 SnpEff
        postHocAnalysis_SnpEff(sampleOutputDir,theSampleName)

        # 3.2 SnpSift
        postHocAnalysis_SnpSift(sampleOutputDir,theSampleName)
        
    return None

def postHocAnalysis_SnpEff(sampleOutputDir,theSampleName):

    '''
    this post hoc analysis for SnpEff.
    '''

    cmd1='%s -Xmx2g -jar %s -v -o gatk -s %s/snps.annotated.summary.html sce.BY4741.toronto %s/%s_specific_snps.vcf > %s/%s_specific_snps.annotated.vcf'%(javaPath,snpEffPath,sampleOutputDir,tmpDir,theSampleName,tmpDir,theSampleName)

    cmd2='%s -Xmx2g -jar %s -v -o gatk -s %s/indels.annotated.summary.html sce.BY4741.toronto %s/%s_specific_indels.vcf > %s/%s_specific_indels.annotated.vcf'%(javaPath,snpEffPath,sampleOutputDir,tmpDir,theSampleName,tmpDir,theSampleName)

    print cmd1
    print
    os.system(cmd1)
    print
    os.system(cmd2)

    return None

def postHocAnalysis_SnpSift(sampleOutputDir,theSampleName):

    path2script=snpEffPath.split('/snpEff.jar')[0]

    cmd1='cat %s/%s_specific_snps.annotated.vcf | perl %s/scripts/vcfEffOnePerLine.pl | %s -jar %s/SnpSift.jar extractFields - CHROM POS REF ALT AF AC DP MQ "(FILTER = \'PASS\')" "ANN[*].ALLELE" "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].GENE" "ANN[*].GENEID" "ANN[*].FEATURE" "ANN[*].BIOTYPE" "ANN[*].RANK" "ANN[*].HGVS_C" "ANN[*].HGVS_P" "ANN[*].CDNA_POS" "ANN[*].CDNA_LEN" "ANN[*].CDS_POS" "ANN[*].CDS_LEN" "ANN[*].AA_POS" "ANN[*].AA_LEN" "ANN[*].DISTANCE" "ANN[*].ERRORS"  > %s/%s/SNPs_finalProduct.txt'%(tmpDir,theSampleName,path2script,javaPath,path2script,outputDir,theSampleName)

    cmd2='cat %s/%s_specific_indels.annotated.vcf | perl %s/scripts/vcfEffOnePerLine.pl | %s -jar %s/SnpSift.jar extractFields - CHROM POS REF ALT AF AC DP MQ "(FILTER = \'PASS\')" "ANN[*].ALLELE" "ANN[*].EFFECT" "ANN[*].IMPACT" "ANN[*].GENE" "ANN[*].GENEID" "ANN[*].FEATURE" "ANN[*].BIOTYPE" "ANN[*].RANK" "ANN[*].HGVS_C" "ANN[*].HGVS_P" "ANN[*].CDNA_POS" "ANN[*].CDNA_LEN" "ANN[*].CDS_POS" "ANN[*].CDS_LEN" "ANN[*].AA_POS" "ANN[*].AA_LEN" "ANN[*].DISTANCE" "ANN[*].ERRORS"  > %s/%s/indels_finalProduct.txt'%(tmpDir,theSampleName,path2script,javaPath,path2script,outputDir,theSampleName)

    print cmd1
    print
    os.system(cmd1)
    print
    print cmd2
    os.system(cmd2)
    print

    return None
        

def sequentialVariantCaller(element):

    pathID=element.split('_1')[0].split('/')[-1]

    tmpGenomePath=tmpDir+'/'+genomePath.split('/')[-1]
    tmpFirstPairFile=tmpDir+'/'+pathID+'/'+element.split('/')[-1]
    tmpSecondPairFile=tmpDir+'/'+pathID+'/mates.fastq'
    workingFile=tmpDir+'/'+pathID+'/file'
    tmpVariantsPath=tmpDir+'/'+variantsPath.split('/')[-1]

    # 1. mapping reads. run bwa with -M (Mark shorter splits as secondary) option in 4 threads (-t)
    print
    print '\t mapping reads...'
    print
    brokenID=element.split('/')[-1].split('_')
    ID='_'.join(brokenID[:3])
    lane=brokenID[3]
    sample=brokenID[4]

    cmd="%s mem -M -t %s -R '@RG\\tID:%s\\tPL:Illumina\\tLB:%s\\tSM:%s' %s %s %s > %s.sam"%(bwaPath,numberOfThreads,ID,lane,sample,tmpGenomePath,tmpFirstPairFile,tmpSecondPairFile,workingFile)
    os.system(cmd)

    # 2. converting SAM to BAM
    print
    print '\tconverting SAM into BAM...'
    print

    cmd='%s -Xmx2g -jar %s SamFormatConverter I=%s.sam O=%s.bam'%(javaPath,piccardPath,workingFile,workingFile)
    os.system(cmd)

    # 3. sorting BAM file
    print
    print '\tsorting BAM file...'
    print
    cmd='%s -Xmx2g -jar %s SortSam I=%s.bam O=%s_sorted.bam SO=coordinate'%(javaPath,piccardPath,workingFile,workingFile)
    os.system(cmd)

    # 4. marking duplicates
    print
    print '\tmarking duplicates...'
    print
    cmd='%s -Xmx2g -jar %s MarkDuplicates I=%s_sorted.bam O=%s_dedup.bam M=%s_duplicatesMetricsInfo.txt'%(javaPath,piccardPath,workingFile,workingFile,workingFile)
    os.system(cmd)

    # 5. building index for bam file
    print
    print '\tbuilding index for bam file...'
    print
    cmd='%s -Xmx2g -jar %s BuildBamIndex I=%s_dedup.bam'%(javaPath,piccardPath,workingFile)
    os.system(cmd)

    # 6. realigning around the indels
    print
    print '\trealigning around the indels...'
    print
    cmd1='%s -Xmx2g -jar %s -nct %s -fixMisencodedQuals -T RealignerTargetCreator -R %s -I %s_dedup.bam -o %s_forIndelRealigner.intervals'%(javaPath,gatkPath,numberOfThreads,tmpGenomePath,workingFile,workingFile)
    cmd2='%s -Xmx2g -jar %s -nct %s -fixMisencodedQuals -T IndelRealigner -R %s -I %s_dedup.bam -known %s -targetIntervals %s_forIndelRealigner.intervals -o %s_realigned.bam'%(javaPath,gatkPath,numberOfThreads,tmpGenomePath,workingFile,tmpVariantsPath,workingFile,workingFile)

    print '\t\tfirst step...'
    print
    os.system(cmd1)
    print
    print '\t\tsecond step...'
    print
    os.system(cmd2)

    # 8. recalibrating BAM file
    print
    print '\tBAM recalibration...'
    print
    
    cmd1='%s -Xmx2g -jar %s -nct %s -T BaseRecalibrator -R %s -I %s_realigned.bam -knownSites %s -o %s_recal_data.grp'%(javaPath,gatkPath,numberOfThreads,tmpGenomePath,workingFile,tmpVariantsPath,workingFile)
    cmd2='%s -Xmx2g -jar %s -nct %s -T PrintReads -R %s -I %s_realigned.bam -BQSR %s_recal_data.grp -o %s_recalibrated.bam'%(javaPath,gatkPath,numberOfThreads,tmpGenomePath,workingFile,workingFile,workingFile)

    print '\t\tfirst step...'
    print
    os.system(cmd1)
    print
    print '\t\tsecond step...'
    print
    os.system(cmd2)

    # 9. sequential variant calling
    print
    print '\tsequential variant calling...'
    print
    cmd='%s -Xmx2g -jar %s -nct %s -T HaplotypeCaller -R %s -I %s_recalibrated.bam -ERC GVCF --variant_index_type LINEAR --variant_index_parameter 128000 -o %s_output.raw.snps.indels.g.vcf'%(javaPath,gatkPath,numberOfThreads,tmpGenomePath,workingFile,workingFile)
    os.system(cmd)
    
    return None

# 0. initialization, basically defining variables and copying files
print 'Welcome to the GATK-bases variant calling pipeline...'

# 0.1 user defined variables
javaPath='/usr/bin/java'
gatkPath='/Users/alomana/software/gatk/GenomeAnalysisTK-3.4-0/GenomeAnalysisTK.jar'
bwaPath='/Users/alomana/software/bwa-0.7.12/bwa'
piccardPath='/Users/alomana/software/picard-tools-1.135/picard.jar'
samtoolsPath='/usr/local/bin/samtools'
snpEffPath='/Users/alomana/software/snpEff/snpEff.jar'

dataPath='/Users/alomana/projects/ap/seqs/data' # this should be a directory with all the read files you want to analyse
genomePath='/Users/alomana/gDrive2/projects/ap/src/gatk/genome/BY4741_Toronto_2012.fasta' # path to the fa file of the genome reference file
variantsPath='/Users/alomana/gDrive2/projects/ap/src/gatk/genome/BY4741_SNV_INDEL_reordered.vcf'

outputDir='/Volumes/WINDOW/ap/gatk'+'/output'
tmpDir='/Volumes/WINDOW/ap/gatk'+'/tmp'

numberOfThreads=4

# 0.2. establishing the directory structure
print
print '\t initializing directory structure...'
print
firstPairs=directoryStructureInitializer()

# 1. analysis pipeline
pipeline(firstPairs)

# 2. post hoc analysis
print
print '\t proceeding into post hoc analysis...'
print
sampleSpecificAnalysis()

print '... done.'
 
