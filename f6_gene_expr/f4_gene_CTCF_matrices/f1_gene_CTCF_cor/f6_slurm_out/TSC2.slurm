#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TSC2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TSC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TSC22D1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TSC22D2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TSC22D3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TSC22D4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TSEN15 
