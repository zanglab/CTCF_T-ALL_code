#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RPL31.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RPL31 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPL31P11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPL32 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPL32P3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPL34 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RPL35 
