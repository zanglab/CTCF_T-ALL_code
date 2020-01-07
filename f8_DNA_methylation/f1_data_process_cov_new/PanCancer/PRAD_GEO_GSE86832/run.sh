#tail -n +2 GSE86832_bigTable.tsv |awk '{OFS="\t";print$1,$2,$2+1,$3,$4,$5,$6}' > GSE86832_bigTable.bed
#liftOver GSE86832_bigTable_hg19.bed /nv/vol190/zanglab/zw5j/data/liftover/hg19ToHg38.over.chain GSE86832_bigTable_hg38.bed unmapped

cat GSE86832_bigTable_hg38.bed|awk '{OFS="\t";if($5>5) print$1,$2,$3,100*$4/$5}' > GSE86832_bigTable_hg38_LNCaP_cov5.bdg
bedtools sort -i GSE86832_bigTable_hg38_LNCaP_cov5.bdg > GSE86832_bigTable_hg38_LNCaP_cov5_sorted.bdg

cat GSE86832_bigTable_hg38.bed|awk '{OFS="\t";if($7>5) print$1,$2,$3,100*$6/$7}' > GSE86832_bigTable_hg38_PrEC_cov5.bdg
bedtools sort -i GSE86832_bigTable_hg38_PrEC_cov5.bdg > GSE86832_bigTable_hg38_PrEC_cov5_sorted.bdg

