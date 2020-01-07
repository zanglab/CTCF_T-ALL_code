#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o IRF3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i IRF3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IRF4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IRF5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IRF6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IRF7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IRF8 
