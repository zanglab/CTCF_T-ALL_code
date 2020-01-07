#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MXI1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MXI1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MXRA5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MXRA7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MXRA8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MYADM 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MYADML 
