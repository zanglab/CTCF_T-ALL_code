#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ANP32C.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ANP32C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANP32D 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANP32E 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANPEP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANTXR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ANTXR2 
