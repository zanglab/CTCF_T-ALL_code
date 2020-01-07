#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MMACHC.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MMACHC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MMADHC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MMD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MMD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MME 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MMEL1 
