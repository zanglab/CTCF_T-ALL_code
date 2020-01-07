#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KIF5A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KIF5A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KIF5B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KIF5C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KIF6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KIF7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KIF9 
