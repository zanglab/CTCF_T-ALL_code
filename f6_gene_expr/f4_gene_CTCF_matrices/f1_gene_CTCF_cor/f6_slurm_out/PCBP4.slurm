#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PCBP4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PCBP4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCCA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCCB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCDH1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCDH10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCDH11X 
