#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o IGHMBP2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i IGHMBP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IGIP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IGLL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IGLL3P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IGLL5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IGLON5 
