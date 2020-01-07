#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NEIL2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NEIL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NEIL3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NEK1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NEK10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NEK11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NEK2 
