#!/bin/bash
#SBATCH -n 1
#SBATCH --mem=80000
#SBATCH -t 24:00:00
#SBATCH -p largemem
#SBATCH -A zanglab 
#SBATCH -o run4.out


ctcf_bam_file="/nv/vol190/zanglab/zw5j/projects_data/T_ALL_CTCF_panos/CTCF/TALL_bam_merge/Jurkat_merge.bam"
atac_bed_file="Jurkat_ATAC_PEnoM.bed"
region_file_mid="union_CTCF_No_Jurkat_DNAmethylation_with_motif_ATAC_peak_overlapped_mid_position.bed"
region_file="union_CTCF_No_Jurkat_DNAmethylation_with_motif_ATAC_peak_overlapped.bed"

python get_RPKM_on_regions_readcount.py -i $region_file_mid -t $ctcf_bam_file -s hg38 -f bam -m -e 200 -o f4_patterns/ctcf_RPKM_e200.csv
python get_RPKM_on_regions_readcount.py -i $region_file_mid -t $atac_bed_file -s hg38 -f bed -m -e 200 -o f4_patterns/atac_RPKM_e200.csv

python get_RPKM_on_regions_readcount.py -i $region_file -t $ctcf_bam_file -s hg38 -f bam -o f4_patterns/ctcf_RPKM.csv
python get_RPKM_on_regions_readcount.py -i $region_file -t $atac_bed_file -s hg38 -f bed -o f4_patterns/atac_RPKM.csv

# python get_pattern_near_site_readcount_write_out.py -s hg38 -f bam -m -w 1000 -b 200 -i ${region_file_mid} -t ${ctcf_bam_file} -o f4_patterns/ctcf_binding_pattern_w1000.csv
# python get_pattern_near_site_readcount_write_out.py -s hg38 -f bed -m -w 1000 -b 200 -i ${region_file_mid} -t ${atac_bed_file} -o f4_patterns/atac_binding_pattern_w1000.csv
