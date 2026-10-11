#!/bin/bash
#SBATCH --job-name=bwa_index
#SBATCH --time=2:00:00
#SBATCH --mem=8gb
#SBATCH --cpus-per-task=1
#SBATCH --output=bwa_index.out
#SBATCH --error=bwa_index.err

# load in bwa module and dependencies
module load gcc/12.1.0 
module load bwa/0.7.17

# move to directory where reference genome files are located
cd /gpfs/data/rmoran-lab/projects/parental_care_gabrielle/Ecrag_ref_NCBI/

bwa index Ecragini.fasta
