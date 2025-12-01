#!/bin/bash
#SBATCH --job-name=use_cutadapt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=48:00:00
#SBATCH --mem=24gb
#SBATCH --output=cutadapt.out
#SBATCH --error=cutadapt.err

#   Set the folder where the analysis is being done
export HOME=/gpfs/data/rmoran-lab/projects/parental_care_gabrielle
mkdir -p "${HOME}/CutAdapt_Out"

#   Specify the adapters used
#   Library prep completed using PerkinElmer NEXTFLEX Rapid XP DNA-Seq Kit HT and Unique Dual Index Barcodes
#   Sequenced using: Illumina NovaSeq 6000 S4 - 2x150 v1.5
#   Barcoded Adapters manual: https://perkinelmer-appliedgenomics.com/wp-content/uploads/2021/11/514150-NEXTflex-Unique-Dual-Index-Barcodes-21-06-SetA_1121.pdf
#   Library Prep Kit manual: https://perkinelmer-appliedgenomics.com/wp-content/uploads/marketing/NEXTFLEX/rapid_xp/NEXTFLEX-Rapid-XP-v2-DNA-Seq-Kit-AG072205_24_MA-V23.08_PE.pdf
#   These are sequences that match the Truseq Dual Index Library (https://teichlab.github.io/scg_lib_structs/methods_html/Illumina.html)
export R1_adapter="GATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNNNATCTCGTATGCCGTCTTCTGCTTG"
export R2_adapter="AATGATACGGCGACCACCGAGATCTACACNNNNNNNNACACTCTTTCCCTACACGACGCTCTTC"



#Load the software needed (parallel, cutadapt, and dependencies)
module load gcc/12.1.0
module load parallel/20240622
module load cutadapt/4.5

#   Define a bash function (a small script within a script) for doing the processing
#   We want to do this because we need to do the same thing for each file
#   Later, we basically tell the computer to do this script for every file
#   Our function is called parallelcut
parallelcut() {
    #Get sample ID from SampleID (first argument supplied to parallel with -a)
    SAMP=${1}

    #Use cutadapt to cut the adapters off
    #-e tells it the maximum error rate
    #-a, -g tells it the adapter to cut (we set the adapter above) (-a is 3' adapter -g is 5' adapter)
    #The input file goes at the end of the command
    #Specify the output file after the > for standard output stream of running a script
    cutadapt -e 0.01 --match-read-wildcards -a ${R1_adapter} -o ${HOME}/CutAdapt_Out/${SAMP}_adtrim_R1.fastq.gz ${HOME}/${SAMP}_R1.fastq.gz
    cutadapt -e 0.01 --match-read-wildcards -a ${R2_adapter} -o ${HOME}/CutAdapt_Out/${SAMP}_adtrim_R2.fastq.gz ${HOME}/${SAMP}_R2.fastq.gz
} 


#   Export function so we can call it with parallel
export -f parallelcut

#Navigate to the directory where we can access our reads and output folder
cd ${HOME}

#Run our bash function to use cut adapt on all files in SampleID.txt
#Use the R1_Adapter.txt and R2_Adapter.txt in the function as well
parallel -j ${SLURM_CPUS_PER_TASK:-8} parallelcut :::: "${HOME}/SampleID.txt"
