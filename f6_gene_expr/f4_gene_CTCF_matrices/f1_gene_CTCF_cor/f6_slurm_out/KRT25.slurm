#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KRT25.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KRT25 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRT26 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRT27 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRT28 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRT3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRT31 
