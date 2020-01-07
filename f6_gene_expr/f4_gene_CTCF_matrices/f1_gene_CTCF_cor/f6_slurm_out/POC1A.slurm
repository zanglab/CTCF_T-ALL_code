#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o POC1A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i POC1A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POC1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POC1B-GALNT4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POC5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PODN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PODNL1 
