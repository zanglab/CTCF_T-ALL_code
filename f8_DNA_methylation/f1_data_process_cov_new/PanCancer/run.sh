# cp /nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/12_new_gained_lost_from_all_features/f5_DNA_methylation/PanCancer/f0_ENCODE_data/unionbedg/*cov5* ENCODE_hg38_sorted/


# time python DNAme_cov2bgd.py
# cp /nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/12_new_gained_lost_from_all_features/f5_DNA_methylation/PanCancer/f0_NCBI_data/process_bdg/HCT116_bdg/GSM1465024_ResultCount_C166TACXX_1_WIT1251A53.hg19_rCRSchrm.fa.realign.mdups.recal.cpg.filtered.sort.CG.6plus2_cov5.bed HCT116_bdg_hg19/

#######
# convert the HCT116 files
#######

# infile="/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/12_new_gained_lost_from_all_features/f5_DNA_methylation/PanCancer/f0_NCBI_data/data_HTC116/GSM1465024_ResultCount_C166TACXX_1_WIT1251A53.hg19_rCRSchrm.fa.realign.mdups.recal.cpg.filtered.sort.CG.6plus2.bed"
# 
# for j in 5 10 
# do 
# 	base_name=$(basename $infile .bed)
# 	echo $base_name $j
# 	cat $infile |awk '{OFS="\t";if($8>'$j') print $1,$2,$3,$7}'>HCT116_bdg/${base_name}_cov${j}.bed
# done




#######
## liftOver
#######

infile="GSM1465024_ResultCount_C166TACXX_1_WIT1251A53.hg19_rCRSchrm.fa.realign.mdups.recal.cpg.filtered.sort.CG.6plus2_cov5.bed"
liftOver HCT116_bdg_hg19/$infile /nv/vol190/zanglab/zw5j/data/liftover/hg19ToHg38.over.chain HCT116_hg38/${infile}_hg38.bed unmapped

infile="GSM2579642_HCT116.control.WGBS.run1_val_1.fq_bismark_bt2_pe.bismark.cov_countfilter_5.bdg"
liftOver HCT116_bdg_hg19/$infile /nv/vol190/zanglab/zw5j/data/liftover/hg19ToHg38.over.chain HCT116_hg38/${infile}_hg38.bed unmapped

infile="GSM257964x_HCT116.control.WGBS.bdg"
liftOver HCT116_bdg_hg19/$infile /nv/vol190/zanglab/zw5j/data/liftover/hg19ToHg38.over.chain HCT116_hg38/${infile}_hg38.bed unmapped



### sort

infile="GSM1465024_ResultCount_C166TACXX_1_WIT1251A53.hg19_rCRSchrm.fa.realign.mdups.recal.cpg.filtered.sort.CG.6plus2_cov5.bed_hg38.bed"
bedtools sort -i HCT116_hg38/$infile > HCT116_hg38_sorted/${infile}_sorted.bdg

infile="GSM2579642_HCT116.control.WGBS.run1_val_1.fq_bismark_bt2_pe.bismark.cov_countfilter_5.bdg_hg38.bed"
bedtools sort -i HCT116_hg38/$infile > HCT116_hg38_sorted/${infile}_sorted.bdg

infile="GSM257964x_HCT116.control.WGBS.bdg_hg38.bed"
bedtools sort -i HCT116_hg38/$infile > HCT116_hg38_sorted/${infile}_sorted.bdg




