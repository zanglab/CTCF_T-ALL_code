#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CSAG3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CSAG3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CSAG4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CSDC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CSDE1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CSE1L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CSF1 
