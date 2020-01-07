#!/bin/bash
#SBATCH -n 8
#SBATCH --mem=20000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o slurm_out/erythroid_progenitor_cells_GSM1816076_H3K27ac.out 

module load bwa
module load samtools

bwa aln  -q 5 -l 32 -k 2 -t 8 /scratch/zw5j/data/index/bwa_index/hg38/hg38.fa /nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f9_histone_modification/download//erythroid_progenitor_cells_GSM1816076_H3K27ac.fastq.gz > /nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f9_histone_modification/download/data_process/process_qc_out/erythroid_progenitor_cells_GSM1816076_H3K27ac/erythroid_progenitor_cells_GSM1816076_H3K27ac_erythroid_progenitor_cells_GSM1816076_H3K27ac.sai
bwa samse /scratch/zw5j/data/index/bwa_index/hg38/hg38.fa /nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f9_histone_modification/download/data_process/process_qc_out/erythroid_progenitor_cells_GSM1816076_H3K27ac/erythroid_progenitor_cells_GSM1816076_H3K27ac_erythroid_progenitor_cells_GSM1816076_H3K27ac.sai /nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f9_histone_modification/download//erythroid_progenitor_cells_GSM1816076_H3K27ac.fastq.gz > /nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f9_histone_modification/download/data_process/process_qc_out/erythroid_progenitor_cells_GSM1816076_H3K27ac/erythroid_progenitor_cells_GSM1816076_H3K27ac_erythroid_progenitor_cells_GSM1816076_H3K27ac.sam
samtools view -bS -q 1 -@ 8 /nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f9_histone_modification/download/data_process/process_qc_out/erythroid_progenitor_cells_GSM1816076_H3K27ac/erythroid_progenitor_cells_GSM1816076_H3K27ac_erythroid_progenitor_cells_GSM1816076_H3K27ac.sam > /nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f9_histone_modification/download/data_process/process_qc_out/erythroid_progenitor_cells_GSM1816076_H3K27ac/erythroid_progenitor_cells_GSM1816076_H3K27ac_erythroid_progenitor_cells_GSM1816076_H3K27ac.bam

macs2 callpeak --SPMR -B -q 0.01 --keep-dup 1 -g hs -t \
/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f9_histone_modification/download/data_process/process_qc_out/erythroid_progenitor_cells_GSM1816076_H3K27ac/erythroid_progenitor_cells_GSM1816076_H3K27ac_erythroid_progenitor_cells_GSM1816076_H3K27ac.bam -n /nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f9_histone_modification/download/data_process/process_qc_out/erythroid_progenitor_cells_GSM1816076_H3K27ac/erythroid_progenitor_cells_GSM1816076_H3K27ac

bdg2bw /nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f9_histone_modification/download/data_process/process_qc_out/erythroid_progenitor_cells_GSM1816076_H3K27ac/erythroid_progenitor_cells_GSM1816076_H3K27ac_treat_pileup.bdg /nv/vol190/zanglab/zw5j/data/Genome/hg38/hg38.chrom.sizes

python 2_chipseq_qc_se.py -i name_list_chipseq_se.txt -o slurm_out/erythroid_progenitor_cells_GSM1816076_H3K27ac_qc.log -s hg38 -b erythroid_progenitor_cells_GSM1816076_H3K27ac -f .fastq.gz