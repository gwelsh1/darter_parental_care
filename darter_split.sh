#!/bin/bash
#SBATCH --job-name=split_reads
#SBATCH --array=1-98
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=12:00:00
#SBATCH --mem=10gb
#SBATCH --output=logs/split_%A_%a.out
#SBATCH --error=logs/split_%A_%a.err

module load gcc/11.3.0
module load samtools/1.18

export HOME=/gpfs/data/rmoran-lab/projects/parental_care_gabrielle

MAP_DIR=${HOME}/Aligned_Reads_Processed/mapped
UNMAP_DIR=${HOME}/Aligned_Reads_Processed/unmapped

mkdir -p ${MAP_DIR}
mkdir -p ${UNMAP_DIR}

SAMP=$(sed -n "${SLURM_ARRAY_TASK_ID}p" SampleID.txt)

echo "Processing ${SAMP} at $(date)"

INPUT=${HOME}/Aligned_Reads_Processed/${SAMP}_RG_markdup.bam 

if [ ! -f ${INPUT} ]; then
    echo "ERROR: Input file not found: ${INPUT}"
    exit 1
fi

echo "Splitting reads for ${SAMP}..."

samtools view -hb -@ 8 -f 4 ${INPUT} > ${UNMAP_DIR}/${SAMP}_unmapped.bam

samtools view -hb -@ 8 -F 260 ${INPUT} > ${MAP_DIR}/${SAMP}_mapped.bam

samtools index ${MAP_DIR}/${SAMP}_mapped.bam

echo "Finished ${SAMP} at $(date)"







