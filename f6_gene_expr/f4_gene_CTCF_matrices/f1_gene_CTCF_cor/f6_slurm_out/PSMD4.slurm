#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PSMD4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PSMD4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PSMD5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PSMD6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PSMD7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PSMD8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PSMD9 
