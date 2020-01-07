#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SERPINB8.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SERPINB8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SERPINB9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SERPINC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SERPIND1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SERPINE1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SERPINE2 
