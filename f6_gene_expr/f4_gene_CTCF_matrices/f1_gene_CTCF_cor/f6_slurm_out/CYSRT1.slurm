#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CYSRT1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CYSRT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CYSTM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CYTH1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CYTH2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CYTH3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CYTH4 
