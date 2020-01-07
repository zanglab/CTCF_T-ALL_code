#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RSC1A1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RSC1A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RSF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RSG1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RSL1D1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RSL24D1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RSPH1 
