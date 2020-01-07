#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o HMOX1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i HMOX1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HMOX2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HMSD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HMX1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HMX2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HMX3 
