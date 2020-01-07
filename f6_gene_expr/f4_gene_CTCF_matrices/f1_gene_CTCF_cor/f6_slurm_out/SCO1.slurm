#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SCO1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SCO1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SCO2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SCOC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SCP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SCP2D1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SCPEP1 
