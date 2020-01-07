#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o BRSK2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i BRSK2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BRWD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BRWD3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BSCL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BSDC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BSG 
