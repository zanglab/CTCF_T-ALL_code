#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CDYL.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CDYL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDYL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CEACAM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CEACAM16 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CEACAM18 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CEACAM19 
