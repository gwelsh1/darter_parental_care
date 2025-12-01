#!/bin/bash
#SBATCH --job-name=align_bwa
#SBATCH --array=1-98%10
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --time=12:00:00
#SBATCH --mem=32gb
#SBATCH --output=logs/bwa_%A_%a.out
#SBATCH --error=logs/bwa_%A_%a.err

module load GCC/13.2.0 
module load seqtk/1.4
module load BWA/0.7.18
module load SAMtools/1.21

export HOME=/gpfs/data/rmoran-lab/projects/parental_care_gabrielle
export REF=/gpfs/data/rmoran-lab/projects/parental_care_gabrielle/Ecrag_ref_NCBI/Ecragini.fasta  # Update this path

mkdir -p ${HOME}/Aligned_Reads

SAMP=$(sed -n "${SLURM_ARRAY_TASK_ID}p" SampleID.txt)

echo "Processing ${SAMP}"

# Interleave reads, align with BWA, convert to sorted BAM
seqtk mergepe \
    ${HOME}/Trimmomatic_Out/${SAMP}_adtrim_trim_pair_R1.fastq.gz \
    ${HOME}/Trimmomatic_Out/${SAMP}_adtrim_trim_pair_R2.fastq.gz \
    | bwa mem -t 16 -k 12 -M -p ${REF} - \
    | samtools view -hbu - \
    | samtools sort -o ${HOME}/Aligned_Reads/${SAMP}_sorted.bam -

echo "Finished ${SAMP}"

