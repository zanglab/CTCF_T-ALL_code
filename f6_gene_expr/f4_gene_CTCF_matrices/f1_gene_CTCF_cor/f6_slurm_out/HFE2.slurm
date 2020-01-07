#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o HFE2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i HFE2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HFM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HGD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HGF 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HGFAC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HGH1 
