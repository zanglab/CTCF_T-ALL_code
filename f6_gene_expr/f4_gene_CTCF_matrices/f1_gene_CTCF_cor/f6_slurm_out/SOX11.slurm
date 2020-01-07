#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SOX11.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SOX11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SOX12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SOX13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SOX14 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SOX15 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SOX17 
