#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PSG1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PSG1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PSG10P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PSG11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PSG2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PSG3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PSG4 
