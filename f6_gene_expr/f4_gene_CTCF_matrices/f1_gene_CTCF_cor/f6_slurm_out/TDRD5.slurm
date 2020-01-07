#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TDRD5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TDRD5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TDRD6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TDRD7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TDRD9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TDRKH 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TDRP 
