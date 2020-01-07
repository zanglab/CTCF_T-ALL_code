#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o HYPK.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i HYPK 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HYPM 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IAH1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IAPP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IARS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IARS2 
