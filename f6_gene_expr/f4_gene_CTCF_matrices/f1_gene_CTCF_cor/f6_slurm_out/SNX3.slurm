#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SNX3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SNX3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SNX30 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SNX31 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SNX32 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SNX33 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SNX4 
