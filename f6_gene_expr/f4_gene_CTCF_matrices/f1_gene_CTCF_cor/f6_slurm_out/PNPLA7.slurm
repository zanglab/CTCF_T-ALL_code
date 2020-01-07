#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PNPLA7.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PNPLA7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PNPLA8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PNPO 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PNPT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PNRC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PNRC2 
