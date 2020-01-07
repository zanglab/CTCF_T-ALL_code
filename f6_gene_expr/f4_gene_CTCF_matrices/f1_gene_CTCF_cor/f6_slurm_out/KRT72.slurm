#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KRT72.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KRT72 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRT73 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRT74 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRT75 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRT76 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRT77 
