#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KRT13.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KRT13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRT14 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRT15 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRT16 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRT16P1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRT16P2 
