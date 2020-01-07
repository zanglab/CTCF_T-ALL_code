#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KIF25.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KIF25 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KIF26A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KIF26B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KIF27 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KIF2A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KIF2B 
