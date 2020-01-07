#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o IDO2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i IDO2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IDS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IDUA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IER2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IER3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IER3IP1 
