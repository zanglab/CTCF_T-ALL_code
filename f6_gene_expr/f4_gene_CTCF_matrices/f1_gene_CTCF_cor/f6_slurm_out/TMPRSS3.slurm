#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TMPRSS3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TMPRSS3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMPRSS4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMPRSS5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMPRSS6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMPRSS7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMPRSS9 
