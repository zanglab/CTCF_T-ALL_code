#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GLS2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GLS2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GLT1D1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GLT6D1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GLT8D1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GLT8D2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GLTP 
