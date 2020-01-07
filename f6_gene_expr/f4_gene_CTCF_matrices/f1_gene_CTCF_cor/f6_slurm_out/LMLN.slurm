#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o LMLN.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i LMLN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LMNA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LMNB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LMNB2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LMNTD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LMNTD2 
