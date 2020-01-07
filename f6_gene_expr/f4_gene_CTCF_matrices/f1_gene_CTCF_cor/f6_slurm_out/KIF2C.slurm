#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KIF2C.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KIF2C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KIF3A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KIF3B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KIF3C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KIF4A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KIF4B 
