#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=40000
#SBATCH -t 24:00:00
#SBATCH -p gpu
#SBATCH -A zanglab
#SBATCH -o slurm_patient.out


pyfile="get_pattern_near_site_readcount_write_out.py"
union_bed="union_binding_occupancy_score_GT3_mid_position.bed"
patient_bam_dir="/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/7_panCancer_DNAmethylation/f2_callpeak/f1_ctcf_callpeak_chilin"

mkdir f2_CTCF_pattern_csv_patient
time python $pyfile -s hg38 -f bam -m -w 1000 -b 200 -i $union_bed -t ${patient_bam_dir}/PD9_CTCF1/attic/PD9_CTCF1_treatment.bam -o f2_CTCF_pattern_csv_patient/PD9_CTCF1.txt
time python $pyfile -s hg38 -f bam -m -w 1000 -b 200 -i $union_bed -t ${patient_bam_dir}/PD9_CTCF2/attic/PD9_CTCF2_treatment.bam -o f2_CTCF_pattern_csv_patient/PD9_CTCF2.txt
time python $pyfile -s hg38 -f bam -m -w 1000 -b 200 -i $union_bed -t ${patient_bam_dir}/PD31_CTCF1/attic/PD31_CTCF1_treatment.bam -o f2_CTCF_pattern_csv_patient/PD31_CTCF1.txt
time python $pyfile -s hg38 -f bam -m -w 1000 -b 200 -i $union_bed -t ${patient_bam_dir}/PD31_CTCF2/attic/PD31_CTCF2_treatment.bam -o f2_CTCF_pattern_csv_patient/PD31_CTCF2.txt
