#union_file="/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f2_union_CTCFs/fz_union_combination/f2_tmp/union_binding_tmp.csv"
#tail -n +2 $union_file |awk -F',' '{OFS="\t";print$1,$NF,$NF+1,$4}' > union_binding_mid_position.bed


python find_overlap_keep_info_NOT_sep_strand.py -a union_binding_mid_position.bed -b hg38_exons.bed -s hg38 -p union_exons_overlapped.bed
python find_overlap_keep_info_NOT_sep_strand.py -a union_binding_mid_position.bed -b hg38_introns.bed -s hg38 -p union_introns_overlapped.bed
python find_overlap_keep_info_NOT_sep_strand.py -a union_binding_mid_position.bed -b hg38_4k_promoter_geneID.bed -s hg38 -p union_promoter_overlapped.bed
