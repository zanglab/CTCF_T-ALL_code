#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SERPINB2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SERPINB2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SERPINB3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SERPINB4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SERPINB5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SERPINB6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SERPINB7 
