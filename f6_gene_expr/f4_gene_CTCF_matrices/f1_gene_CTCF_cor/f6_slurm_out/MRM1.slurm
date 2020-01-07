#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MRM1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MRM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MRO 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MROH1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MROH2A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MROH2B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MROH5 
