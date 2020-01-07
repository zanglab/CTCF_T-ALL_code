#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o EDN2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i EDN2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EDN3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EDNRA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EDNRB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EDRF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EEA1 
