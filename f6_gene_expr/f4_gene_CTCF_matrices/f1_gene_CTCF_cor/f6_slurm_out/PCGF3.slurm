#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PCGF3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PCGF3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCGF5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCGF6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCID2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCIF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PCK1 
