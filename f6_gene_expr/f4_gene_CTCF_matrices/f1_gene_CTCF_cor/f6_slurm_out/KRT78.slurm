#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KRT78.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KRT78 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRT79 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRT8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRT80 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRT81 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRT82 
