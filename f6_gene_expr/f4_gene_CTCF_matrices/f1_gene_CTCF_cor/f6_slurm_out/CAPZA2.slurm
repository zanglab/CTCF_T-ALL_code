#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CAPZA2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CAPZA2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CAPZA3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CAPZB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CARD10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CARD11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CARD14 
