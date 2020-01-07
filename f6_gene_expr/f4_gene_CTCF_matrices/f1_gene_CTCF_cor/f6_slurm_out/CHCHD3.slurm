#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CHCHD3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CHCHD3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHCHD4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHCHD5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHCHD6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHCHD7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CHD1 
