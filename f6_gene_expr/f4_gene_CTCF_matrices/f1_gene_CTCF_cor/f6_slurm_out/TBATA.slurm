#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TBATA.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TBATA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TBC1D1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TBC1D10A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TBC1D10B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TBC1D10C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TBC1D12 
