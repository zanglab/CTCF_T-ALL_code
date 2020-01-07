
hic_boundary_merged="f0_HiC_TADs_boundaries/hic_domain_boundary_merge_200k_hg38.bed"

union_bed="/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f2_union_CTCFs/fz_union_combination/f2_union_binding_data/union_binding_occupancy_score_GT3.bed"
consitutive_bed="/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/updated_201906/f2_union_CTCFs/fz_union_combination/f2_union_binding_data/union_binding_constitutive_0.8.bed"
cousitutive_boundary="all_CTCF_domainInfo_raw_boundaries.bed"

# get merged bed file
python bedfile_merge_revised.py -a $cousitutive_boundary -p all_CTCF_domainInfo_raw_boundaries_merged.bed -s hg38 -e 40000
python bedfile_merge_revised.py -a $consitutive_bed -p constitutive_merged.bed -s hg38 -e 40000
python bedfile_merge_revised.py -a $union_bed -p union_merged.bed -s hg38 -e 40000

cousitutive_boundary_merged="all_CTCF_domainInfo_raw_boundaries_merged.bed"
consitutive_merged="constitutive_merged.bed"
union_merged="union_merged.bed"


python find_overlap_keep_info_NOT_sep_strand.py -s hg38 -a $union_bed -b $hic_boundary_merged -p overlapped/union_overlapped_hic_boundaries.bed 
python find_overlap_keep_info_NOT_sep_strand.py -s hg38 -a $union_merged -b $hic_boundary_merged -p overlapped/union_merged_overlapped_hic_boundaries.bed 
python find_overlap_keep_info_NOT_sep_strand.py -s hg38 -a $consitutive_bed -b $hic_boundary_merged -p overlapped/constitutive_overlapped_hic_boundaries.bed 
python find_overlap_keep_info_NOT_sep_strand.py -s hg38 -a $consitutive_merged -b $hic_boundary_merged -p overlapped/constitutive_merged_overlapped_hic_boundaries.bed 
python find_overlap_keep_info_NOT_sep_strand.py -s hg38 -a $cousitutive_boundary -b $hic_boundary_merged -p overlapped/const_boundary_overlapped_hic_boundaries.bed 
python find_overlap_keep_info_NOT_sep_strand.py -s hg38 -a $cousitutive_boundary_merged -b $hic_boundary_merged -p overlapped/const_boundary_merged_overlapped_hic_boundaries.bed 

python find_overlap_keep_info_NOT_sep_strand.py -s hg38 -b $union_bed -a $hic_boundary_merged -p overlapped/hic_boundaries_overlapped_union.bed 
python find_overlap_keep_info_NOT_sep_strand.py -s hg38 -b $union_merged -a $hic_boundary_merged -p overlapped/hic_boundaries_overlapped_union_merged.bed 
python find_overlap_keep_info_NOT_sep_strand.py -s hg38 -b $consitutive_bed -a $hic_boundary_merged -p overlapped/hic_boundaries_overlapped_constitutive.bed 
python find_overlap_keep_info_NOT_sep_strand.py -s hg38 -b $consitutive_merged -a $hic_boundary_merged -p overlapped/hic_boundaries_overlapped_constitutive_merged.bed 
python find_overlap_keep_info_NOT_sep_strand.py -s hg38 -b $cousitutive_boundary -a $hic_boundary_merged -p overlapped/hic_boundaries_overlapped_const_boundary.bed 
python find_overlap_keep_info_NOT_sep_strand.py -s hg38 -b $cousitutive_boundary_merged -a $hic_boundary_merged -p overlapped/hic_boundaries_overlapped_const_boundary_merged.bed 
