#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DCTN6.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DCTN6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DCTPP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DCUN1D1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DCUN1D2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DCUN1D3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DCUN1D4 
