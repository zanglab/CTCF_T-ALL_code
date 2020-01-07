#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o HAPLN1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i HAPLN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HAPLN2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HAPLN3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HAPLN4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HARBI1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HARS 
