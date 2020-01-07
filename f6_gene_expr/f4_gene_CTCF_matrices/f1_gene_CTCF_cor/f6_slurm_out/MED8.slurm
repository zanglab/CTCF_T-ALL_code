#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MED8.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MED8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MED9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MEDAG 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MEF2A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MEF2B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MEF2C 
