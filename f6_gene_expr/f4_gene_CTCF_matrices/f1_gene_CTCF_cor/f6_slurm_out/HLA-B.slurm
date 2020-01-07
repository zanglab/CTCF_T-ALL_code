#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o HLA-B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i HLA-B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HLA-C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HLA-DMA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HLA-DMB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HLA-DOA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HLA-DOB 
