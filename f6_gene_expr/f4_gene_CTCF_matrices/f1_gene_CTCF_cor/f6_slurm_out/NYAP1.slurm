#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NYAP1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NYAP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NYAP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NYNRIN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NYX 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OACYLP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OAF 
