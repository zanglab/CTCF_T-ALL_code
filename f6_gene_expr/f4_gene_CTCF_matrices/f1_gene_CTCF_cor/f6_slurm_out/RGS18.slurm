#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RGS18.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RGS18 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RGS19 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RGS2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RGS20 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RGS21 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RGS22 
