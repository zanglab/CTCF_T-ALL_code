#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CABP1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CABP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CABP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CABP4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CABP5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CABP7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CABS1 
