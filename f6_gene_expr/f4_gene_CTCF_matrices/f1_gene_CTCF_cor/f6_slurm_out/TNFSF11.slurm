#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TNFSF11.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TNFSF11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TNFSF12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TNFSF12-TNFSF13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TNFSF13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TNFSF13B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TNFSF14 
