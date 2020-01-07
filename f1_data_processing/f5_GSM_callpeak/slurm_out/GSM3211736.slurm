#!/bin/bash
#SBATCH -n 4
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o slurm_out/GSM3211736.out 

module load bwa
module load samtools

bwa aln  -q 5 -l 32 -k 2 -t 8 /scratch/zw5j/data/index/bwa_index/hg38/hg38.fa /nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f1_data_processing/f3_fq_files/GSM3211736.fastq.gz > /nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f1_data_processing/f5_GSM_callpeak/GSM_callpeak_out/GSM3211736/GSM3211736.sai
bwa samse /scratch/zw5j/data/index/bwa_index/hg38/hg38.fa /nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f1_data_processing/f5_GSM_callpeak/GSM_callpeak_out/GSM3211736/GSM3211736.sai /nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f1_data_processing/f3_fq_files/GSM3211736.fastq.gz > /nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f1_data_processing/f5_GSM_callpeak/GSM_callpeak_out/GSM3211736/GSM3211736.sam
samtools view -bS -q 1 -@ 8 /nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f1_data_processing/f5_GSM_callpeak/GSM_callpeak_out/GSM3211736/GSM3211736.sam > /nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f1_data_processing/f5_GSM_callpeak/GSM_callpeak_out/GSM3211736/GSM3211736.bam

macs2 callpeak --SPMR -B -q 0.01 --keep-dup 1 -g hs -t \
/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f1_data_processing/f5_GSM_callpeak/GSM_callpeak_out/GSM3211736/GSM3211736.bam -n /nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f1_data_processing/f5_GSM_callpeak/GSM_callpeak_out/GSM3211736/GSM3211736

bdg2bw /nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f1_data_processing/f5_GSM_callpeak/GSM_callpeak_out/GSM3211736/GSM3211736_treat_pileup.bdg /nv/vol190/zanglab/zw5j/data/Genome/hg38/hg38.chrom.sizes

python 2_chipseq_qc_se_macs2_tags_no_reps.py -i name_list_chipseq_se.txt -o slurm_out/GSM3211736_qc.log -s hg38 -b GSM3211736 -f .fastq.gz