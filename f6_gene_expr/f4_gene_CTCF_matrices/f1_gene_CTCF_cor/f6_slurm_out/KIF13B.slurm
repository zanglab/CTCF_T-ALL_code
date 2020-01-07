#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KIF13B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KIF13B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KIF14 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KIF15 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KIF16B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KIF17 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KIF18A 
