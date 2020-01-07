#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CTAGE15.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CTAGE15 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CTAGE4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CTAGE5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CTAGE6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CTAGE7P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CTAGE8 
