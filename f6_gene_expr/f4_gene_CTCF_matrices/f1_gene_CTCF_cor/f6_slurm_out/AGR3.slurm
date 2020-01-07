#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o AGR3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i AGR3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AGRN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AGRP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AGT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AGTPBP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AGTR1 
