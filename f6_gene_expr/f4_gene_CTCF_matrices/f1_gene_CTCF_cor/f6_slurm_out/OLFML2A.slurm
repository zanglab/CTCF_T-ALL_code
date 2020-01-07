#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o OLFML2A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i OLFML2A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OLFML2B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OLFML3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OLIG1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OLIG2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OLIG3 
