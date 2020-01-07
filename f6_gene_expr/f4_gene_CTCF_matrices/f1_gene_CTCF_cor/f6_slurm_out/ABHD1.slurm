#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ABHD1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ABHD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ABHD10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ABHD11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ABHD11-AS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ABHD12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ABHD12B 
