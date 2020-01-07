#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ADCY10P1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ADCY10P1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADCY2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADCY3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADCY4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADCY5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADCY6 
