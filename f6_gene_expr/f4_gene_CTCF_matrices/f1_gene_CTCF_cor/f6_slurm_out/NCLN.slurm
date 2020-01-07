#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NCLN.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NCLN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NCMAP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NCOA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NCOA2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NCOA3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NCOA4 
