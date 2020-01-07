#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DDX11L16.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DDX11L16 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DDX11L2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DDX11L5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DDX11L9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DDX12P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DDX17 
