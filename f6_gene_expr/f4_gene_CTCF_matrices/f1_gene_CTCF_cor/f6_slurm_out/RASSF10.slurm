#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RASSF10.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RASSF10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RASSF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RASSF3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RASSF4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RASSF5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RASSF6 
