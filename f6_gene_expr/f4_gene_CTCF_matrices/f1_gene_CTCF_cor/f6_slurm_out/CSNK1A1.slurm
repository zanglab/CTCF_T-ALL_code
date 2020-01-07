#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CSNK1A1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CSNK1A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CSNK1A1L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CSNK1A1P1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CSNK1D 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CSNK1E 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CSNK1G1 
