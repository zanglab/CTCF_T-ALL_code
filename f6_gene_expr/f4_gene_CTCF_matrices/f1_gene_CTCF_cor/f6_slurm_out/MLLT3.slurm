#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MLLT3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MLLT3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MLLT6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MLN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MLNR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MLPH 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MLST8 
