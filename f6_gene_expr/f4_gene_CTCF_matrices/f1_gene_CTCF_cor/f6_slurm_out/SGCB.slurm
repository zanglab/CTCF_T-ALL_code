#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SGCB.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SGCB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SGCD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SGCE 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SGCG 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SGCZ 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SGF29 
