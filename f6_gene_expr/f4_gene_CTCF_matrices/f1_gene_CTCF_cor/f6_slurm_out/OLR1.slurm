#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o OLR1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i OLR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OMA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OMD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OMG 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OMP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ONECUT1 
