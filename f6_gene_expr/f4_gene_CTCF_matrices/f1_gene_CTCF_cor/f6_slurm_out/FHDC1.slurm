#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FHDC1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FHDC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FHIT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FHL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FHL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FHL3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FHL5 
