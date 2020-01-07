#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MTG2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MTG2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MTHFD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MTHFD1L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MTHFD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MTHFD2L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MTHFD2P1 
