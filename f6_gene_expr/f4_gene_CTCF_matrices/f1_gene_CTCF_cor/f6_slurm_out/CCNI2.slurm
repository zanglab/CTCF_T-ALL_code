#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CCNI2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CCNI2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCNJ 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCNJL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCNK 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCNL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CCNL2 
