#!/bin/bash
#SBATCH --job-name=index_bams
#SBATCH --array=1-98
#SBATCH --cpus-per-task=1
#SBATCH --mem=4gb
#SBATCH --time=1:00:00
#SBATCH --output=logs/index_%A_%a.out
#SBATCH --error=logs/index_%A_%a.err

module load gcc/11.3.0
module load samtools/1.18

export HOME=/gpfs/data/rmoran-lab/projects/parental_care_gabrielle

# Get sample name
SAMP=$(sed -n "${SLURM_ARRAY_TASK_ID}p" SampleID.txt)

echo "Indexing ${SAMP}"

samtools index ${HOME}/Aligned_Reads/${SAMP}_sorted.bam

echo "Done indexing ${SAMP}"
ls -lh ${HOME}/Aligned_Reads/${SAMP}_sorted.bam*
EOF