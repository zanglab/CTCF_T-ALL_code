#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TBC1D3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TBC1D3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TBC1D30 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TBC1D31 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TBC1D32 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TBC1D3B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TBC1D3C 
