#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PTH2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PTH2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTH2R 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTHLH 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTK2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTK2B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PTK6 
