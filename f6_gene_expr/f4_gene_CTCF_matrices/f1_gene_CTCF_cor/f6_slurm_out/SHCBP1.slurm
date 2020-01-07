#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SHCBP1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SHCBP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SHCBP1L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SHD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SHE 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SHF 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SHH 
