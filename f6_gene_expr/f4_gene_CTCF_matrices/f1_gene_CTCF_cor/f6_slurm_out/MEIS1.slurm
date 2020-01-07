#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MEIS1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MEIS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MEIS2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MEIS3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MEIS3P1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MELK 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MEMO1 
