#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o VARS2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i VARS2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VASH1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VASH2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VASN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VASP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VAT1 
