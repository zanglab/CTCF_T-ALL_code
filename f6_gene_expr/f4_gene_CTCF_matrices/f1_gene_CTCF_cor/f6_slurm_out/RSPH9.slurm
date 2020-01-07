#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RSPH9.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RSPH9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RSPO1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RSPO2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RSPO3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RSPO4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RSPRY1 
