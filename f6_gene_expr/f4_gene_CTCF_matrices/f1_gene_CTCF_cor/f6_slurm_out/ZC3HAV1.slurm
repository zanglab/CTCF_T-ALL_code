#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZC3HAV1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZC3HAV1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZC3HAV1L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZC3HC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZC4H2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZCCHC10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZCCHC11 
