# Activating Qiime2
source activate qiime2-2018.11

# Adding paths and directories
DATADIR="fastq_files"
Primers="/Metagenetics_Projects/Primers/Primers_16SV3-V4_18S-V4.txt" # This file contains the 16S and 18S forward and reverse primers. Make sure the last line of the file is empty to avoid skipping of the last line.
mkdir Demultiplexed_fastq_files
Output="Demultiplexed_fastq_files"
mkdir $Output/fastq_16S
mkdir $Output/fastq_18S

# Loop going through each file and looking for the forward and reverse primers using cutadapt implemented in Qiime2
# Primers are trimmed in the demultiplexed files

for i in $DATADIR/*_R1_001.fastq.gz
do

R1=$i;
R2=${R1%_R1_001.fastq.gz}_R2_001.fastq.gz;
SAMPLE_F=`basename ${R1};`
SAMPLE_R=`basename ${R2};`

while read -r Primer_name Forward_primers Reverse_primers; do
cutadapt -g X$Forward_primers -G X$Reverse_primers -o $Output/$Primer_name$SAMPLE_F -p $Output/$Primer_name$SAMPLE_R $R1 $R2 -j 0 --no-indels --discard-untrimmed --overlap 17

done<$Primers

# Write down the summary of the demultiplexing process
done > demultiplexing_summary.txt

# Export the demultiplexed files into new directories
mv $Output/16S* $Output/fastq_16S
mv $Output/18S* $Output/fastq_18S
