#!/bin/bash
#SBATCH -n 4
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o slurm_out//GSM3212930.out 

module load bwa
module load samtools

bwa aln  -q 5 -l 32 -k 2 -t 8 /scratch/zw5j/data/index/bwa_index/hg38/hg38.fa /nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f1_data_processing/f3_fq_files/GSM3212930_1.fastq.gz > /nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f1_data_processing/f5_GSM_callpeak/GSM_callpeak_out/GSM3212930/GSM3212930.1.sai
bwa aln  -q 5 -l 32 -k 2 -t 8 /scratch/zw5j/data/index/bwa_index/hg38/hg38.fa /nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f1_data_processing/f3_fq_files/GSM3212930_2.fastq.gz > /nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f1_data_processing/f5_GSM_callpeak/GSM_callpeak_out/GSM3212930/GSM3212930.2.sai
bwa sampe /scratch/zw5j/data/index/bwa_index/hg38/hg38.fa /nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f1_data_processing/f5_GSM_callpeak/GSM_callpeak_out/GSM3212930/GSM3212930.1.sai /nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f1_data_processing/f5_GSM_callpeak/GSM_callpeak_out/GSM3212930/GSM3212930.2.sai /nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f1_data_processing/f3_fq_files/GSM3212930_1.fastq.gz /nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f1_data_processing/f3_fq_files/GSM3212930_2.fastq.gz > /nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f1_data_processing/f5_GSM_callpeak/GSM_callpeak_out/GSM3212930/GSM3212930.sam
samtools view -bS -q 1 -@ 8 /nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f1_data_processing/f5_GSM_callpeak/GSM_callpeak_out/GSM3212930/GSM3212930.sam > /nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f1_data_processing/f5_GSM_callpeak/GSM_callpeak_out/GSM3212930/GSM3212930.bam

macs2 callpeak --SPMR -B -q 0.01 --keep-dup 1 -g hs -f BAMPE -t \
/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f1_data_processing/f5_GSM_callpeak/GSM_callpeak_out/GSM3212930/GSM3212930.bam -n /nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f1_data_processing/f5_GSM_callpeak/GSM_callpeak_out/GSM3212930/GSM3212930

bdg2bw /nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f1_data_processing/f5_GSM_callpeak/GSM_callpeak_out/GSM3212930/GSM3212930_treat_pileup.bdg /nv/vol190/zanglab/zw5j/data/Genome/hg38/hg38.chrom.sizes

python 2_chipseq_qc_pe_macs2_fragments_no_reps.py -i name_list_chipseq_pe.txt -o slurm_out//GSM3212930_qc.log -s hg38 -b GSM3212930 -f1 _1.fastq.gz -f2 _2.fastq.gz