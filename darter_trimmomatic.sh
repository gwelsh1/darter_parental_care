#!/bin/bash
#SBATCH --job-name=trimmomatic
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=48:00:00
#SBATCH --mem=24gb
#SBATCH --output=trimmomatic.out
#SBATCH --error=trimmomatic.err

#   Set the folder where the analysis is being done
export HOME=/gpfs/data/rmoran-lab/projects/parental_care_gabrielle
mkdir -p "${HOME}/Trimmomatic_Out"
#Input files in ./CutAdapt_Out
#Output files in ./Trimmomatic_Out

#       Load software modules needed
module load gcc/12.1.0
module load parallel/20240622
module load openjdk/11.0.2
module load openjdk/17.0.2
module load trimmomatic/0.39

#NOTE:

#   Define a bash function (a small script within a script) for doing the processing
#   We want to do this because we need to do the same thing for each file
#   Later, we basically tell the computer to do this script for every file
#   Our function is called paralleltrim
paralleltrim() {
        #Get each sample name from the list SampleID.txt (supplied to the parallel command with the -a flag)
        SAMP=${1}
        
        #trim low quality
        java -Xmx2g -jar /apps/software/openjdk-17.0.2/trimmomatic/0.39/trimmomatic-0.39.jar PE -phred33 ${HOME}/CutAdapt_Out/${SAMP}_adtrim_R1.fastq.gz ${HOME}/CutAdapt_Out/${SAMP}_adtrim_R2.fastq.gz ${HOME}/Trimmomatic_Out/${SAMP}_adtrim_trim_pair_R1.fastq.gz ${HOME}/Trimmomatic_Out/${SAMP}_adtrim_trim_unpair_R1.fastq.gz ${HOME}/Trimmomatic_Out/${SAMP}_adtrim_trim_pair_R2.fastq.gz ${HOME}/Trimmomatic_Out/${SAMP}_adtrim_trim_unpair_R2.fastq.gz SLIDINGWINDOW:6:30 MINLEN:40
}

#   Export function so we can call it with parallel
export -f paralleltrim

#NAvigate into the reads directory
cd ${HOME}

parallel --joblog trimmomatic_parallel_logfile.txt -a SampleID.txt paralleltrim


