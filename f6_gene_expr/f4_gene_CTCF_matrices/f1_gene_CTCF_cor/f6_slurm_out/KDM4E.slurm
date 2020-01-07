#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KDM4E.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KDM4E 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KDM5A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KDM5B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KDM5C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KDM5D 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KDM6A 
