#!/bin/bash
#SBATCH --job-name=make_picardjobs
#SBATCH --array=1-98
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --mem=10gb
#SBATCH --time=2:00:00
#SBATCH --output=logs/picard_%A_%a.out
#SBATCH --error=logs/picard_%A_%a.err

module load openjdk/17.0.2 
module load picard/3.4.0

export HOME=/gpfs/data/rmoran-lab/projects/parental_care_gabrielle

mkdir -p ${HOME}/Aligned_Reads_Processed
mkdir -p logs

SAMP=$(sed -n "${SLURM_ARRAY_TASK_ID}p" SampleID.txt)

echo "Processing ${SAMP} at $(date)"

java -Xmx4G -jar ${PICARD} AddOrReplaceReadGroups \
	TMP_DIR=/tmp \
    INPUT=${HOME}/Aligned_Reads/${SAMP}_sorted.bam \
    OUTPUT=${HOME}/Aligned_Reads_Processed/${SAMP}_RG.bam \
    RGID=${SAMP} \
    RGLB=${SAMP} \
    RGPL=ILLUMINA \
    RGPU=${SAMP} \
    RGSM=${SAMP} \
    CREATE_INDEX=true

echo "Read groups added for ${SAMP}"

java -Djava.io.tmpdir=${TMPDIR} -Xmx12g -jar ${PICARD} MarkDuplicates \
    TMP_DIR=${TMPDIR} \
    INPUT=${HOME}/Aligned_Reads_Processed/${SAMP}_RG.bam \
    OUTPUT=${HOME}/Aligned_Reads_Processed/${SAMP}_RG_markdup.bam \
    METRICS_FILE=${HOME}/Aligned_Reads_Processed/${SAMP}_dup_metrics.txt \
    REMOVE_DUPLICATES=false \
    CREATE_INDEX=true

echo "Finished ${SAMP} at $(date)"

if [ -f ${HOME}/Aligned_Reads_Processed/${SAMP}_dup_metrics.txt ]; then
    echo "Duplicate metrics:"
    tail -n 5 ${HOME}/Aligned_Reads_Processed/${SAMP}_dup_metrics.txt
fi