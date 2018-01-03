#!/bin/bash -l

#assuming availability of
#bowtie2
#SAMtools
#python 3 (for python requirements see requirements.txt)

#set the working directory to the directory in which this script resides
SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
  TARGET="$(readlink "$SOURCE")"
  if [[ $TARGET == /* ]]; then
    echo "SOURCE '$SOURCE' is an absolute symlink to '$TARGET'"
    SOURCE="$TARGET"
  else
    working_dir="$( dirname "$SOURCE" )"
    echo "SOURCE '$SOURCE' is a relative symlink to '$TARGET' (relative to '$working_dir')"
    SOURCE="$working_dir/$TARGET" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
  fi
done
echo "SOURCE is '$SOURCE'"
RDIR="$( dirname "$SOURCE" )"
working_dir="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
if [ "$working_dir" != "$RDIR" ]; then
  echo "working_dir '$RDIR' resolves to '$working_dir'"
fi
echo "working_dir is '$working_dir'"

RAW=${working_dir}/SRA_CAGE_reads.fastq.gz
SAMPLE_BARCODE=${working_dir}/sample_barcode.csv

#check the path to the AGPv3 bowtie2 index
#to get this file and folder download the genome_index.tar.gz file in the working dir
REFERENCE=${working_dir}'/genome_index/AGPv3/ZmB73'

#PREFIX_TAGS holds SRA strings to build files names
PREFIX_TAGS=(\
SRR2078285 \
SRR2078286 \
SRR2078287 \
SRR2078288 \
SRR2078289 \
SRR2078290 \
SRR2078291 \
SRR2078292)

#READ_GROUP_ARRAY holds the id for the read group
SAMPLE_ARRAY=(\
B73_shoot \
B73_shoot \
B73_root \
B73_root \
Mo17_shoot \
Mo17_shoot \
Mo17_root \
Mo17_root)

#SAMPLE_RRAY holds the taxa/tissue/bioreplicate information for naming
READ_GROUP_ARRAY=(\
SRR2078285_GATCAGCAG \
SRR2078286_ACACAGCAG \
SRR2078287_ACTCAGCAG \
SRR2078288_ACGCAGCAG \
SRR2078289_AGACAGCAG \
SRR2078290_ATCCAGCAG \
SRR2078291_ATGCAGCAG \
SRR2078292_CTTCAGCAG)


#1a. If a single *.fastq file with all the reads then using the preprocessing_reads.py scripts makes sense
#re-demultiplex reads and trimming to release the CAGE tag
#minimun pre-processing is ideal, as the TSS position could be clipped and lost
#be aware of this step when interpreting the TSSs shifting as some of the conclusions might be artifacts of the preprocessing!
#python preprocessing_reads.py $RAW $SAMPLE_BARCODE

for((i=0;i<${#SAMPLE_ARRAY[@]};i++));
do
	#set up variables
	SAMPLE = ${working_dir}/${PREFIX_TAGS[i]}"_demultiplexed_trimmed_filtered.fastq"
	SORTED_BAM_1 = ${working_dir}/${PREFIX_TAGS[i]}_sorted.bam
    FILTERED_BAM = ${working_dir}/${PREFIX_TAGS[i]}_sorted_mapq_20.bam

	#1b. assuming that files were downloaded from the SRA and converted individually into fastq files
	#reads already demultiplexed and trimmed to release the CAGE tag
	#minimun pre-processing is ideal, as the TSS position could be clipped and lost!
	#echo '*****************trimming '${SAMPLE_ARRAY[i]}'*****************'
	python barcode_sequence_trimmer.py ${PREFIX_TAGS[i]}

    echo '*****************BOWTIE2 '${PREFIX_TAGS[i]}'*****************'
    #2. alignments bowtie2
    #alignment using bowtie2, sorting BAM file with SAMTOOLS (*sorted.bam)
    bowtie2 -p 8 --rg-id ${READ_GROUP_ARRAY[i]} --rg "SM:"${SAMPLE_ARRAY[i]} --rg "LB:CAGE_Grotewold_lab" --rg "PL:ILLUMINA"  --very-sensitive-local -x $REFERENCE $SAMPLE | samtools view -Shu - | samtools sort -O BAM - > $SORTED_BAM_1
    
    #3. indexing alignment output BAM files 
	echo '*****************samtools index alignment output '${PREFIX_TAGS[i]}'*****************'
	#indexing sorted BAM file with SAMTOOLS (*.bai)
    samtools index -b $SORTED_BAM_1
    
    #4. Post-alignment filtering for MAPQ >= 20 BAM file
    echo '*****************samtools filtering sorted BAM '${PREFIX_TAGS[i]}'*****************'
    samtools view -F 4 -q 20 -u $SORTED_BAM_1 | samtools sort -O BAM -@ 10 - > $FILTERED_BAM
    
    #5. indexing filtered BAM files 
	echo '*****************samtools index filtered output '${PREFIX_TAGS[i]}'*****************'
	#indexing sorted BAM file with SAMTOOLS (*.bai)
    samtools index -b $FILTERED_BAM
    
done

