#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o VHLL.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i VHLL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VIL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VILL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VIM 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VIP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VIPAS39 
