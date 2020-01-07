#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PDE11A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PDE11A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PDE12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PDE1A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PDE1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PDE1C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PDE2A 
