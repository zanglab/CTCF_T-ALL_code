#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TBC1D24.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TBC1D24 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TBC1D25 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TBC1D26 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TBC1D28 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TBC1D29 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TBC1D2B 
