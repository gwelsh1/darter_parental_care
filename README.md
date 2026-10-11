# darter_parental_care
This is a repository for scripts used by Gabrielle Welsh to process darter WGS data. Primary goals of the project were identifying areas of the genome under selection in darter lineages that repeatedly gained parental care and under relaxed selection in lineages that secondarily lost that trait.

Step 1: Raw reads processed with CutAdapt and Trimmomatic.
(see 01_cutadapt.sh and 02_trimmomatic.sh)

Step 2: Reads were aligned to a Etheostoma cragini reference genome (downloaded from NCBI)
(see 03_index_ref.sh and 04_align.sh)


