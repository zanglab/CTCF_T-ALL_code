#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o AGPAT2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i AGPAT2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AGPAT3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AGPAT4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AGPAT5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AGPS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AGR2 
