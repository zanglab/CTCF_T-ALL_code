#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MB.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MB21D2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MBD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MBD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MBD3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MBD3L1 
