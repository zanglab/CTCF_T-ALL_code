#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SIRPB1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SIRPB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SIRPB2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SIRPD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SIRPG 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SIRT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SIRT2 
