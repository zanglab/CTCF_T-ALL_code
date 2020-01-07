#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MYDGF.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MYDGF 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MYEF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MYEOV 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MYF5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MYF6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MYH1 
