#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SPANXN3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SPANXN3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPANXN4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPANXN5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPARC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPARCL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPAST 
