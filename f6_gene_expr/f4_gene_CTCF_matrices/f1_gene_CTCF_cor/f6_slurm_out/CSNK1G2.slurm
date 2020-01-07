#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CSNK1G2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CSNK1G2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CSNK1G3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CSNK2A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CSNK2A2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CSNK2A3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CSNK2B 
