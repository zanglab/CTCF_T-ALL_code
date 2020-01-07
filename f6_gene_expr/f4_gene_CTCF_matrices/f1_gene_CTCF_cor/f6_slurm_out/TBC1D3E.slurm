#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TBC1D3E.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TBC1D3E 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TBC1D3F 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TBC1D3G 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TBC1D3H 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TBC1D3I 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TBC1D3K 
