#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SERPINA9.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SERPINA9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SERPINB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SERPINB10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SERPINB11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SERPINB12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SERPINB13 
