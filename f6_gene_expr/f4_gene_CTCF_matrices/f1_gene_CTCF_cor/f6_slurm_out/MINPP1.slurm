#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MINPP1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MINPP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MIOS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MIOX 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MIP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MIPEP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MIPEPP3 
