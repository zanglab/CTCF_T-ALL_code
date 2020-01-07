#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o BAZ1A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i BAZ1A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BAZ1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BAZ2A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BAZ2B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BBC3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BBIP1 
