#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CHEK2P2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CHEK2P2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHERP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHFR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHGA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHGB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHI3L1 
