#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NCR1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NCR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NCR2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NCR3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NCR3LG1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NCS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NCSTN 
