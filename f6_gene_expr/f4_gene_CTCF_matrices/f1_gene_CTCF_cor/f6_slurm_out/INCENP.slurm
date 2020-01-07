#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o INCENP.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i INCENP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i INF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ING1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ING2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ING3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ING4 
