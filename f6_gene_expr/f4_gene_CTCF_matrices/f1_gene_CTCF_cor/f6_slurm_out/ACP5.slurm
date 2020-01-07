#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ACP5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ACP5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACP6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACPP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACRBP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACRV1 
