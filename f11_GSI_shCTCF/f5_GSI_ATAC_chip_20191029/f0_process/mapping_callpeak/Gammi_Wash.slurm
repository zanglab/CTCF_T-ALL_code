#!/bin/bash
#SBATCH -n 2
#SBATCH --mem=60000
#SBATCH -t 24:00:00
#SBATCH -p largemem
#SBATCH -A zanglab
#SBATCH -o Gammi_Wash.out


Trim5=15
Trim3=5
S="hg38"
macs2_S="hs"

fq_dir="/nv/vol190/zanglab/shared/TALL_mayo/Transfer/"
prename="C1_Gammi_Wash"

mkdir ${prename}
O=${prename}/${prename}

bowtie2 -p 8 -5 ${Trim5} -3 ${Trim3} -X 2000 -x /nv/vol190/zanglab/sh8tv/Data/Mapping_index/bowtie2/${S}_clean -1 ${fq_dir}/${prename}_R1.fastq.gz -2 ${fq_dir}/${prename}_R2.fastq.gz -S ${O}.sam 2>&1 >>/dev/null | tee ${O}_bowtie2PE.out

samtools view -bS ${O}.sam > ${O}_raw.bam
# filtering & samTObam
samtools view -H ${O}_raw.bam > ${O}_filtered.sam 
samtools view -f 0x2 ${O}_raw.bam | awk 'NR % 2 == 1{mapq=$5;forward=$0} NR % 2 == 0{if($5>=30 && mapq>=30 && substr($3,1,3)=="chr" && $7=="=" ) print forward"\n"$0}' >> ${O}_filtered.sam 
samtools view -bS ${O}_filtered.sam > ${O}.bam
wait
rm ${O}.sam ${O}_filtered.sam 

# PE process
bamToBed -i ${O}.bam -bedpe | awk '{OFS="\t"; print $1,$2,$6,$7,$8,"."}' | sort -k 1,1 -k 2,2g -k 3,3g > ${O}_PE.bed 
grep chrM ${O}_PE.bed > ${O}_PEchrM.bed
grep -v chrM ${O}_PE.bed > ${O}_PEnoM.bed

macs2 callpeak --SPMR -B -q 0.01 --keep-dup 1 -f BAMPE -g ${macs2_S} -t ${O}.bam -n ${O}
mv ${O}_treat_pileup.bdg ${O}.bdg
bdg2bw ${O}.bdg /nv/vol190/zanglab/sh8tv/Data/Genome/${S}/${S}_clean.len



