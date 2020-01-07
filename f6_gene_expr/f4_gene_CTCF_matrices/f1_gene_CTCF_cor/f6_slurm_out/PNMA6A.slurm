#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PNMA6A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PNMA6A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PNMT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PNN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PNO1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PNOC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PNP 
