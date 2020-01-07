#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MSS51.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MSS51 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MST1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MST1L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MST1P2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MST1R 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MSTN 
