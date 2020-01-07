
for i in ../domain_deg/*gained_up*GT100K_LT1M_log2FC0.585_padj1e-3.txt
do 
  bart2 geneset -i ${i} -s hg38 --outdir bart2_result
done

for i in ../domain_deg/*lost_down*GT100K_LT1M_log2FC0.585_padj1e-3.txt
do 
  bart2 geneset -i ${i} -s hg38 --outdir bart2_result
done
