#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CCNC.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CCNC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCND1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCND2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCND3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCNDBP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCNE1 
