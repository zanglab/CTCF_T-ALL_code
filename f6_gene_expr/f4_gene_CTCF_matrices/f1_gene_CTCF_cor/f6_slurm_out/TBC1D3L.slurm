#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TBC1D3L.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TBC1D3L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TBC1D3P2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TBC1D3P5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TBC1D4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TBC1D5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TBC1D7 
