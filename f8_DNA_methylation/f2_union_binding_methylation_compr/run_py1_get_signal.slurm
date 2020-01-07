#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=40000
#SBATCH -t 24:00:00
#SBATCH -p gpu
#SBATCH -A zanglab 
#SBATCH -o run1.out


# ==== PanCancer ====

motif_bed_file="/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f2_union_CTCFs/fz_union_combination/f2_union_binding_data/union_binding_occupancy_score_GT3_mid_position.bed"
pancancer_methylation_dir="/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f8_DNA_methylation/f1_data_process_cov_new/PanCancer/ENCODE_hg38_sorted"
    
time python get_loci_signal_on_regions_scores.py -i $motif_bed_file -t ${pancancer_methylation_dir}/A549_cov5_unioned.bdg  -o f1_signal_binding_csv/A549_cov5_DNAme_unionBinding.csv -s hg38 -c 0 -m -e 150
time python get_loci_signal_on_regions_scores.py -i $motif_bed_file -t ${pancancer_methylation_dir}/lung_cov5_unioned.bdg  -o f1_signal_binding_csv/lung_cov5_DNAme_unionBinding.csv -s hg38 -c 0 -m -e 150
time python get_loci_signal_on_regions_scores.py -i $motif_bed_file -t ${pancancer_methylation_dir}/lung_cov5_unioned_WGBS.bdg  -o f1_signal_binding_csv/lung_cov5_WGBS_DNAme_unionBinding.csv -s hg38 -c 0 -m -e 150

time python get_loci_signal_on_regions_scores.py -i $motif_bed_file -t ${pancancer_methylation_dir}/MCF7_cov5_unioned.bdg  -o f1_signal_binding_csv/MCF7_cov5_DNAme_unionBinding.csv -s hg38 -c 0 -m -e 150
time python get_loci_signal_on_regions_scores.py -i $motif_bed_file -t ${pancancer_methylation_dir}/breast_cov5_unioned.bdg  -o f1_signal_binding_csv/breast_cov5_DNAme_unionBinding.csv -s hg38 -c 0 -m -e 150

time python get_loci_signal_on_regions_scores.py -i $motif_bed_file -t ${pancancer_methylation_dir}/HCT116_cov5_unioned.bdg  -o f1_signal_binding_csv/HCT116_cov5_DNAme_unionBinding.csv -s hg38 -c 0 -m -e 150
time python get_loci_signal_on_regions_scores.py -i $motif_bed_file -t ${pancancer_methylation_dir}/sigmoid_colon_cov5_unioned.bdg  -o f1_signal_binding_csv/sigmoid_colon_cov5_DNAme_unionBinding.csv -s hg38 -c 0 -m -e 150


hct116_file_dir="/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f8_DNA_methylation/f1_data_process_cov_new/PanCancer/HCT116_hg38_sorted"
time python get_loci_signal_on_regions_scores.py -i $motif_bed_file -t ${hct116_file_dir}/GSM2579642_HCT116.control.WGBS.run1_val_1.fq_bismark_bt2_pe.bismark.cov_countfilter_5.bdg_hg38.bed_sorted.bdg  -o f1_signal_binding_csv/HCT116_GSM1465024_cov5_WGBS_DNAme_unionBinding.csv -s hg38 -c 0 -m -e 150
time python get_loci_signal_on_regions_scores.py -i $motif_bed_file -t ${hct116_file_dir}/GSM257964x_HCT116.control.WGBS.bdg_hg38.bed_sorted.bdg  -o f1_signal_binding_csv/HCT116_GSM257964x_cov5_WGBS_DNAme_unionBinding.csv -s hg38 -c 0 -m -e 150



# ==== T-ALL ====

tall_methylation_dir="/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/8_omics_analysis_cancer_specific_CTCF_lost/f2_heatmap_on_omics_data_T_ALL_CTCFlost/f0_infiles/TALL_omics_files/DNAmethylation"

time python get_loci_signal_on_regions_scores.py -i $motif_bed_file -t ${tall_methylation_dir}/JURKAT_DMSO_hg38_union.bdg  -o f1_signal_binding_csv/Jurkat_DNAme_unionBinding.csv -s hg38 -c 0 -m -e 150
time python get_loci_signal_on_regions_scores.py -i $motif_bed_file -t ${tall_methylation_dir}/CUTLL1_hg38_union.bdg  -o f1_signal_binding_csv/CUTLL1_DNAme_unionBinding.csv -s hg38 -c 0 -m -e 150
time python get_loci_signal_on_regions_scores.py -i $motif_bed_file -t ${tall_methylation_dir}/A6_CD4_hg38_union.bdg  -o f1_signal_binding_csv/A6010_DNAme_unionBinding.csv -s hg38 -c 0 -m -e 150



# ==== T-ALL patients ====

patient_methylation_dir="/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f8_DNA_methylation/f1_data_process_cov_new/patient/patient_bdg"
    
time python get_loci_signal_on_regions_scores.py -i $motif_bed_file -t ${patient_methylation_dir}/Pt9_countfilter_5.bdg   -o f1_signal_binding_csv/PD9_cov5_DNAme_unionBinding.csv -s hg38 -c 0 -m -e 150
time python get_loci_signal_on_regions_scores.py -i $motif_bed_file -t ${patient_methylation_dir}/Pt31_countfilter_5.bdg  -o f1_signal_binding_csv/PD31_cov5_DNAme_unionBinding.csv -s hg38 -c 0 -m -e 150




