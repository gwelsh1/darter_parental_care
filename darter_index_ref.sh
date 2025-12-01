#!/bin/bash
#SBATCH --job-name=bwa_index
#SBATCH --time=2:00:00
#SBATCH --mem=8gb
#SBATCH --cpus-per-task=1
#SBATCH --output=bwa_index.out
#SBATCH --error=bwa_index.err

module load gcc/12.1.0 
module load bwa/0.7.17

cd /gpfs/data/rmoran-lab/projects/parental_care_gabrielle/Ecrag_ref_NCBI/


echo "Starting BWA indexing at $(date)..."
bwa index Ecragini.fasta

echo "Indexing complete at $(date)!"
echo "Files created:"
ls -lh Ecragini.fasta.*
EOF