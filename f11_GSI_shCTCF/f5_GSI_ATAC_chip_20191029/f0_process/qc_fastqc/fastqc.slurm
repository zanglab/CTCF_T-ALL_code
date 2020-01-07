#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=30000
#SBATCH -t 24:00:00
#SBATCH -p gpu
#SBATCH -A zanglab
#SBATCH -o qc.out



# module load cutadapt 
# module load fastqc
# export PATH="/nv/vol190/zanglab/zw5j/env/trim_galore/TrimGalore-0.5.0:$PATH"

mkdir fastqc_out

fastqc -o fastqc_out/ /nv/vol190/zanglab/shared/TALL_mayo/Transfer/A1_DMSO_3d_R1.fastq.gz   
fastqc -o fastqc_out/ /nv/vol190/zanglab/shared/TALL_mayo/Transfer/A1_DMSO_3d_R2.fastq.gz   
fastqc -o fastqc_out/ /nv/vol190/zanglab/shared/TALL_mayo/Transfer/B1_Gammi_3d_R1.fastq.gz  
fastqc -o fastqc_out/ /nv/vol190/zanglab/shared/TALL_mayo/Transfer/B1_Gammi_3d_R2.fastq.gz
fastqc -o fastqc_out/ /nv/vol190/zanglab/shared/TALL_mayo/Transfer/C1_Gammi_Wash_R1.fastq.gz
fastqc -o fastqc_out/ /nv/vol190/zanglab/shared/TALL_mayo/Transfer/C1_Gammi_Wash_R2.fastq.gz


# trim_galore --nextera --phred33 --fastqc --paired --clip_R1 30 --clip_R2 30 --three_prime_clip_R1 30 --three_prime_clip_R2 30 SRR891277_1.fastq.gz SRR891277_2.fastq.gz -o test2_trim30/
# trim_galore --nextera --phred33 --fastqc --paired SRR5063984_1.fastq.gz SRR5063984_2.fastq.gz -o trimmed_out
# trim_galore --nextera --phred33 --fastqc --paired SRR5063985_1.fastq.gz SRR5063985_2.fastq.gz -o trimmed_out

