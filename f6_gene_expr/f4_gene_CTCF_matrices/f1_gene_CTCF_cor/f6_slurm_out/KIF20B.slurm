#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KIF20B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KIF20B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KIF21A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KIF21B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KIF22 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KIF23 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KIF24 
