#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SH3PXD2A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SH3PXD2A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SH3PXD2B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SH3RF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SH3RF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SH3RF3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SH3TC1 
