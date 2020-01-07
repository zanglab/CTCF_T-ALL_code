#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FHOD1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FHOD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FHOD3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FIBCD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FIBIN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FIBP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FICD 
