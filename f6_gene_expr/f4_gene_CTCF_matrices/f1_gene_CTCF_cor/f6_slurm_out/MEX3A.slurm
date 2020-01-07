#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MEX3A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MEX3A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MEX3B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MEX3C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MEX3D 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MFAP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MFAP2 
