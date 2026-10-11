#!/bin/bash
#SBATCH --job-name=align_bwa
#SBATCH --array=1-118%10
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --time=12:00:00
#SBATCH --mem=32gb
#SBATCH --output=logs/bwa_%A_%a.out
#SBATCH --error=logs/bwa_%A_%a.err


# load in modules
module load gcc/11.3.0
module load seqtk/1.3
module load bwa/0.7.17
module load samtools/1.18

# designate paths for raw reads and the indexed reference genome
export HOME=/gpfs/data/rmoran-lab/projects/parental_care_gabrielle
export REF=/gpfs/data/rmoran-lab/projects/parental_care_gabrielle/Ecrag_ref_NCBI/Ecragini.fasta  

# make an output directory called "Aligned_Reads"
mkdir -p ${HOME}/Aligned_Reads

# use the SampleID.txt list to define each sample
SAMP=$(sed -n "${SLURM_ARRAY_TASK_ID}p" SampleID.txt)

# Interleave reads, align with BWA, convert to sorted BAM
seqtk mergepe \
    ${HOME}/Trimmomatic_Out/${SAMP}_adtrim_trim_pair_R1.fastq.gz \
    ${HOME}/Trimmomatic_Out/${SAMP}_adtrim_trim_pair_R2.fastq.gz \
    | bwa mem -t 16 -k 12 -M -p ${REF} - \
    | samtools view -hbu - \
    | samtools sort -o ${HOME}/Aligned_Reads/${SAMP}_sorted.bam -




