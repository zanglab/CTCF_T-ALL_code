#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RGS11.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RGS11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RGS12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RGS13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RGS14 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RGS16 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RGS17 
