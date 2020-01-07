#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MICAL3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MICAL3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MICALL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MICALL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MICB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MICU1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MICU2 
