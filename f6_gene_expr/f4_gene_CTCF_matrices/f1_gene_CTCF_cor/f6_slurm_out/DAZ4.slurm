#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DAZ4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DAZ4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DAZAP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DAZAP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DAZL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DBF4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DBF4B 
