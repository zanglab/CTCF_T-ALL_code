#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KRT37.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KRT37 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRT38 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRT39 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRT4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRT40 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KRT42P 
