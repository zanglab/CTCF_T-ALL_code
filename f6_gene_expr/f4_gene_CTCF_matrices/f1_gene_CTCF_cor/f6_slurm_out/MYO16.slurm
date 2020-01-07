#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MYO16.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MYO16 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MYO18A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MYO18B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MYO19 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MYO1A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MYO1B 
