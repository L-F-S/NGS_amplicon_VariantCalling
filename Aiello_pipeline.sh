#!/bin/bash

####################################################################################
# @AUTHOR			Lorenzo Fedrico Signorini		           #
# @PROJECT			 EvolvR eGFP variant calling                       #
# @DESCRIPTION		The main shell script to execute the pipeline		   #
# @DATE                                 09/01/2020                                 #
####################################################################################

# 1) Make sure that data is in $BASE_PATH

# 2) Make sure a fasta file with the reference sequence in .fna format is in
# $REF_PATH

# 3) Set the right path for the following programs:
# [fastqc, trim_galore]

# 4) In order to execute this script you need to allow the execution permissions
# Thus, after locating the file, run: chmod +x comano_pipeline.sh

# PATH VARIABLES DEFINITIONS ###################################################
BASE_PATH='/shares/CIBIO-Storage/CM/scratch/tmp_projects/signorini_cas/Aiello_EvolvR'		# base path of the project
REF_PATH=${BASE_PATH}
DATASET=${BASE_PATH}'/00_raw_data'
PREPROCESSING=${BASE_PATH}'/01_trimmed_data'
MAPPING=${BASE_PATH}'/02_mapping'
VC=${BASE_PATH}'/03_variant_calling'

FASTQC_INIT=${DATASET}'/fastqc_report'
FASTQC_TG=${PREPROCESSING}'/trimgalore_fastqc'
TG_REPORTS=${PREPROCESSING}'/trim_galore_reports'
TM_UNPAIRED=${PREPROCESSING}'/trimmomatic_unpaired'
TM_RESULTS_PLOT=${PREPROCESSING}'/trimmomatic_results_plot'
ANNOTATIONS=${ASSEMBLY}'/annotations'

# TOOLS VARIABLES DEFINITIONS #################################################
trim_galore='/shares/CIBIO-Storage/CM/mir/tools/trim_galore-0.4.2/trim_galore'
TrimmomaticPE='/home/lorenzo.signorini/mapping_CJ_vars/utils/TrimmomaticPE.1'
fastqc='/shares/CIBIO-Storage/CM/mir/tools/fastqc-0.11.8/fastqc'
ref_sequence_file='amplicon_ref'  #avoid extension. make sure is .fna
samtools='/shares/CIBIO-Storage/CM/mir/tools/samtools-1.10/bin/samtools'


echo -e '############################################################\n'
echo -e ' EvolvR eGFP variant calling\n'
echo -e ' Author: L-F-S\n'
echo -e '############################################################\n'
#
#
# INSPECTION OF THE ORIGINAL DATA ##############################################
echo -e '# Running *FastQC* on the original data...\n'
mkdir ${FASTQC_INIT}
$fastqc ${DATASET}/*.fastq.gz -o ${FASTQC_INIT}
echo '# Reports successfully created in:'
echo -e '#     '${FASTQC_INIT}'\n'


# PREPROCESSING OF SEQUENCING DATA #############################################
# This trims lower-quality positions (usually at the end of reads),
# removes adapters (if present),
# and produces fastQC reports of the newly polished data.
echo -e '# Running *trim_galore* on the original data...\n'
mkdir ${PREPROCESSING} && mkdir ${FASTQC_TG} && mkdir ${TG_REPORTS}
#$trim_galore -q 30 --fastqc ${DATASET}/*.fastq.gz -o ${PREPROCESSING}  --path_to_cutadapt /home/lorenzo.signorini/.local/bin/cutadapt
$trim_galore -q 30 -o ${PREPROCESSING}  --path_to_cutadapt /home/lorenzo.signorini/.local/bin/cutadapt
mv ${PREPROCESSING}/*.txt ${TG_REPORTS}
$fastqc ${PREPROCESSING}/*.fq.gz -o ${FASTQC_TG} 
echo '# Data successfully preprocessed by *trim_galore*. Check in:'
echo -e '#     '${PREPROCESSING}'\n'


# DATA MAPPING #########################################
echo -e '# Sequences mapping'
echo -e '# building reference index..'
cd ${REF_PATH}   # to run bowtie, you must be in the directory of the reference
bowtie2-build ${ref_sequence_file}.fna $ref_sequence_file

echo -e '# running *BowTie2* on trimmed data..' 

# UNCOMMENT ONE OF THE TWO:
echo -e '# aligning single-end reads'
for SAMPLE in ${PREPROCESSING}/*_R1_001_trimmed.fq.gz; do  #change output read name
	echo 'sample ' $SAMPLE
	bowtie2 -x ${ref_sequence_file} -U $SAMPLE -S ${SAMPLE}.sam &
        BACK_PID=$!
done

#echo -e '# aligning paired-end reads..'
#for INPUT_FORWARD in ${PREPROCESSING}/*R1_001_*; do
#	echo forward $INPUT_FORWARD
#	INPUT_REVERSE=${INPUT_FORWARD/R1_001_val_1/R2_001_val_2}
#	echo reverse $INPUT_REVERSE
#	N=$(echo ${INPUT_FORWARD} | cut -d '_' -f 9  )
#	echo $N
#	bowtie2 -x $ref_sequence_file -1 ${INPUT_FORWARD} -2 ${INPUT_REVERSE} -S S${N}_CJ_variant.sam &
#	BACK_PID=$!
#done

echo -e 'waiting for alignments to finish'
wait $BACK_PID  #only keeping track of last bowtie of the cycle.l should be 
       		#the last to finish

cd $PREPROCESSING
mkdir $MAPPING
mv *.sam ${MAPPING}
cd ${MAPPING}
echo -e '# Compressing SAM files'
for file in *.sam; do echo $file; tar -czvf ${file}.tar.gz $file & done
echo -e '# Converting SAM to BAM, sorting and indexing...'
cd $MAPPING
for sam in *sam; do 
	echo $sam
	$samtools view -bS $sam > ${sam/sam/bam}  # must use uncompressed sam
	$samtools sort ${sam/sam/bam} -o ${sam/sam/sorted}  # can use compressed\
	                                              #	bam. Sometimes works
#	#	only without -o flag, sometimes only with (?)
	$samtools index ${sam/sam/sorted}
	#$samtools index ${sam/bam.tar.gz/sorted}
done

#echo -e '# Creating per-base depth files'
#for file in *sorted; do
#	depth=${file/sorted/cov}
#	echo -e 'creating file '$depth
#	$samtools depth $file > $depth &
#done

rm *sam
sono stati creati i bam compressi? rm *bam

## VARIANTS CALLING v1 ################################
#echo -e '# Variants calling'
#cd $BASE_PATH
#mkdir 03_vcf
#for N in 0 1 2 3; do
#	SBAM='SS'$N'_CJ_variant.sorted.bam'
#	echo $SBAM
#	OUTVCF='S'$N'_CJ_variant.bcf'
#	samtools mpileup -uf reference/CJ_plasmid.fna 02_mapping/$SBAM | bcftools view -v - > ${VC}/$OUTVCF
#done

#VARIANTS CALLING V2: PACBAM
mkdir $VC
echo -e "# Using PaCBAM"

# creating a bed file
echo -e "# Creating .bed and.vcf files to run paCBAM.."
cd $REF_PATH
SEQ_NAME='Ref_201'
echo -e $SEQ_NAME '\t0\t201'>${ref_sequence_file}.bed
echo -e '# Fake vcf file needed to run PaCBAM'>fake.vcf
echo -e '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT'>>fake.vcf
echo -e $SEQ_NAME '\t1\t.\tG\t.\t.\t.\t.\t.'>>fake.vcf


# not deduped:
cd $MAPPING
for SBAM in *sorted; do
	echo $SBAM
	pacbam bam=$SBAM bed=$REF_PATH/$ref_sequence_file.bed vcf=$REF_PATH/fake.vcf fasta=${REF_PATH}/${ref_sequence_file}.fna out=$VC/ & 
done

 deduped
cd $MAPPING
for SBAM in *sorted; do
	echo $SBAM
	pacbam bam=$SBAM bed=$REF_PATH/$ref_sequence_file.bed vcf=$REF_PATH/fake.vcf fasta=${REF_PATH}/${ref_sequence_file}.fna out=${VC}_dedup/ dedup &
