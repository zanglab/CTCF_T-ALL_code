#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o LIMD2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i LIMD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LIME1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LIMK1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LIMK2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LIMS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LIMS2 
