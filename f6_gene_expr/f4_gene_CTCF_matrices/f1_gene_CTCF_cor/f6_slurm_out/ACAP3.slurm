#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ACAP3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ACAP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACAT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACAT2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACBD3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACBD4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ACBD5 
