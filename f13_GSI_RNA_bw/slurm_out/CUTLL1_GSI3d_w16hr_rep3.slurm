#!/bin/bash
#SBATCH -n 2
#SBATCH --mem=40000
#SBATCH -t 24:00:00
#SBATCH -p gpu
#SBATCH -A zanglab
#SBATCH -o slurm_out/CUTLL1_GSI3d_w16hr_rep3.out 

index_file="/nv/vol190/zanglab/zw5j/data/index/star_index/ucsc_refGene/hg38"
rsem_reference="/nv/vol190/zanglab/zw5j/data/index/rsem_reference/ucsc_refGene/hg38/rsem_reference"
mapping_outdir="processed_results//CUTLL1_GSI3d_w16hr_rep3"

STAR \
--runMode alignReads \
--genomeDir ${index_file} \
--readFilesIn \
/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/12_new_gained_lost_from_all_features/f6_Jurkat_GSI_RNA_seq_new/f0_RNA_trimmed/data_trimmed/trimmed_out//GSI_wo_94_C3_S17_L001_R1_001_trimmed.fq.gz,\
/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/12_new_gained_lost_from_all_features/f6_Jurkat_GSI_RNA_seq_new/f0_RNA_trimmed/data_trimmed/trimmed_out//GSI_wo_94_C3_S17_L002_R1_001_trimmed.fq.gz,\
/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/12_new_gained_lost_from_all_features/f6_Jurkat_GSI_RNA_seq_new/f0_RNA_trimmed/data_trimmed/trimmed_out//GSI_wo_94_C3_S17_L003_R1_001_trimmed.fq.gz,\
/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/12_new_gained_lost_from_all_features/f6_Jurkat_GSI_RNA_seq_new/f0_RNA_trimmed/data_trimmed/trimmed_out//GSI_wo_94_C3_S17_L004_R1_001_trimmed.fq.gz \
 \
--outFileNamePrefix ${mapping_outdir}/CUTLL1_GSI3d_w16hr_rep3 \
--outFilterMultimapNmax 1 \
--runThreadN 4 \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outWigType wiggle \
--quantMode TranscriptomeSAM GeneCounts \
--readFilesCommand zcat

