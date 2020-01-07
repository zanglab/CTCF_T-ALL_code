#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PDK4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PDK4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PDLIM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PDLIM2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PDLIM3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PDLIM4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PDLIM5 
