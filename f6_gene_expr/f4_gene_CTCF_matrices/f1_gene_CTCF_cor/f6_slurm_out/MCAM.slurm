#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MCAM.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MCAM 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MCAT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MCC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MCCC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MCCC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MCCD1 
