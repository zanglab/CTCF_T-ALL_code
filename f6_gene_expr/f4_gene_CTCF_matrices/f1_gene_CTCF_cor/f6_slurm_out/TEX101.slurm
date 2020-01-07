#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TEX101.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TEX101 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TEX11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TEX12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TEX13A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TEX13B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TEX14 
