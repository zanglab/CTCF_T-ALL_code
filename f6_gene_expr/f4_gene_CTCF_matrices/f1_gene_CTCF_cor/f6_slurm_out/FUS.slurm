#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FUS.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FUS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FUT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FUT10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FUT11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FUT2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FUT3 
