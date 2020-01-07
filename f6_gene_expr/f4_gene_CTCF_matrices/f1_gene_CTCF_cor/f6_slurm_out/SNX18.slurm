#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SNX18.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SNX18 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SNX19 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SNX2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SNX20 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SNX21 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SNX22 
