#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PSG5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PSG5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PSG6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PSG7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PSG8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PSG9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PSIP1 
