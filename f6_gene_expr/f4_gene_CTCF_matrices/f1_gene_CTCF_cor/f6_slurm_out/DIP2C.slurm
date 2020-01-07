#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DIP2C.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DIP2C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DIRAS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DIRAS2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DIRAS3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DIRC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DIRC2 
