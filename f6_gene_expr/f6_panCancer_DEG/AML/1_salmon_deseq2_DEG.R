#source("https://bioconductor.org/biocLite.R")
#biocLite("DESeq2")
#biocLite("tximport")
#biocLite("tximportData")
library("tximport")
library("DESeq2")
library("tximportData")
library(gplots)
library(RColorBrewer)
source('plotPCAWithSampleNames.R')
#install.packages("ggrepel")

###################
#COPY REGION START 
################### 

file_dir = getwd()
setwd(file_dir)
tx2gene_dir="/nv/vol190/zanglab/zw5j/data/index/salmon_index/tx2gene_ensembl2symbol"
tx2gene <- read.csv(file.path(tx2gene_dir, "hg38_v92_tx2gene.csv"),header=TRUE,sep="\t")

# read the salmon-out of all samples
salmon_out_dir = "/nv/vol190/zanglab/zw5j/work2017/T_ALL_CTCF/10_CTCF_binding_signals_vs_gene_expression/f4_CTCF_cor_gene_cancer_vs_normal/f1_process_RNA/salmon_quant"

# get the col match info of all samples
sampleTable_tmp = read.table(file.path(file_dir,"colNames.txt"),header=TRUE)
# I do not know why but after rewrite the sampleTable, the result could run correctly
sampleTable <- data.frame(condition=factor(sampleTable_tmp$condition),id=factor(sampleTable_tmp$id),name=factor(sampleTable_tmp$name))
rownames(sampleTable) <- sampleTable_tmp$id

sample_files <- file.path(salmon_out_dir,sampleTable_tmp$id, "quant.sf")
all(file.exists(sample_files))
names(sample_files) <- sampleTable_tmp$id


###################
# COPY REGION END 
###################

txi=tximport(sample_files,type='salmon',tx2gene=tx2gene)
save(txi,file="f1_R_odds_saved/txi")
#load(file="f1_R_odds_saved/txi")
sampleTable <- sampleTable[colnames(txi$counts),]

ddsTxi=DESeqDataSetFromTximport(txi,colData=sampleTable,design=~condition)
dds=DESeq(ddsTxi)
save(dds,file = "f1_R_odds_saved/dds")
#load(file = "f1_R_odds_saved/dds")
#dir.create('f2_samlmon_pca')

log_dds<-rlog(dds)
save(log_dds,file = "f1_R_odds_saved/log_dds")
#load(file = "f1_R_odds_saved/log_dds")


####################
# PACA analysis
####################
library(RColorBrewer)
colors = brewer.pal(9, "Set1")
#colorCodes <- c('1'="red", 'a'="green", C="blue", D="yellow")
labelCol <- function(x) {
  if (is.leaf(x)) {
    ## fetch label
    label <- attr(x, "label")
    code <- sampleTable[label,]$name
    #code <- substr(label, 1, 1);#sprintf(code);sprintf('a');sprintf('a')
    ## use the following line to reset the label to one letter code
    # attr(x, "label") <- code
    attr(x, "nodePar") <- list(lab.col=colors[code])
  }
  return(x)
}


source('plotPCAWithSampleNames.R')
pdf("f2_samlmon_pca/all_samples_pca_log_allgenes.pdf", width=5, height=5)
plotPCAWithSampleNames(log_dds, intgroup="name", ntop=40000,label_col=4)
dev.off()

dists <- dist(t(assay(log_dds)))
hc <- hclust(dists)
d <- dendrapply(as.dendrogram(hc), labelCol)
pdf("f2_samlmon_pca/all_samples_hclust_log_allgenes.pdf", width=6, height=4)
par(mar=c(15,2,0,0)+1)
plot(d)
dev.off()


#####################
# identify DEG
#####################

func_deseq2<-function(sampleTable,txi,dds,treat,ctrl,csv_name)
{
counts_table=txi$abundance
#counts_table <-counts_table[!rowSums(counts_table==0)>=10,]
# keep only those data in a/b group
filter_id = sampleTable[(sampleTable$condition==treat)|(sampleTable$condition==ctrl),]
counts_table = counts_table[,filter_id$id]

filtered_norm_counts<-counts_table[rowSums(counts_table!=0)>=1,]
filtered_norm_counts<-as.data.frame(filtered_norm_counts)
GeneID<-rownames(filtered_norm_counts)
filtered_norm_counts<-cbind(filtered_norm_counts,GeneID)
dim(filtered_norm_counts)
#head(filtered_norm_counts)

#res_shrink <- lfcShrink(dds, coef=paste0("condition_",treat,"_vs_",ctrl), type="apeglm")
res_shrink <- lfcShrink(dds, contrast=c("condition",treat,ctrl))
res_ordered<-res_shrink[order(res_shrink$padj),]
GeneID<-rownames(res_ordered)
res_ordered<-as.data.frame(res_ordered)
res_genes<-cbind(res_ordered,GeneID)
dim(res_genes)
res_genes_merged <- merge(filtered_norm_counts,res_genes,by=unique("GeneID"))
res_ordered<-res_genes_merged[order(res_genes_merged$padj),]
dim(res_ordered)
#head(res_ordered)
write.csv(res_ordered, file = file.path("f4_deseq_out_shrink_saved/",csv_name),row.names=FALSE)
}


CtrlList = c('normal')
TreatList = c('cancer')

for (i in 1:1){
treat = TreatList[i]
ctrl = CtrlList[i]
names_sep = c("treated_",treat,"_vs_ctrl_",ctrl,".csv")
func_deseq2(sampleTable,txi,dds,treat,ctrl,paste(names_sep,collapse=""))
}

#####record the gene level abundance/TMP
abundance_table=txi$abundance
filtered_norm_abundance<-abundance_table[rowSums(abundance_table!=0)>=1,]
filtered_norm_abundance<-as.data.frame(filtered_norm_abundance)
dim(filtered_norm_abundance)
head(filtered_norm_abundance)
write.csv(filtered_norm_abundance, file = file.path("f4_deseq_out_shrink/","gene_level_abundance_TPM_all.csv"),row.names=TRUE)

