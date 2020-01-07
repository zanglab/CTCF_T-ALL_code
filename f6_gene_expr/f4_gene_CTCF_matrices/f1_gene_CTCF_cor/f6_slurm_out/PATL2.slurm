#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PATL2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PATL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PATZ1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PAWR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PAX1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PAX2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PAX3 
