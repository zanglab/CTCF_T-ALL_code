#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ERAP1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ERAP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ERAP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ERAS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ERBB2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ERBB3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ERBB4 
