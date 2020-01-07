#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o EMCN.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i EMCN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EMD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EME1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EME2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EMG1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EMID1 
