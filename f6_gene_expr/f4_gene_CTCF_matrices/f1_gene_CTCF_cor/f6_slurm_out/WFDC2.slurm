#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o WFDC2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i WFDC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WFDC21P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WFDC3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WFDC5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WFDC6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WFDC8 
