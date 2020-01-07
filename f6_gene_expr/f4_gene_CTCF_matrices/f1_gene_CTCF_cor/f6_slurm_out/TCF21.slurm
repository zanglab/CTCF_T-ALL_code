#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TCF21.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TCF21 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TCF23 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TCF24 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TCF25 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TCF3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TCF4 
