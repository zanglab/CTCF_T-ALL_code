#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ENC1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ENC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ENDOD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ENDOG 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ENDOU 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ENDOV 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ENG 
