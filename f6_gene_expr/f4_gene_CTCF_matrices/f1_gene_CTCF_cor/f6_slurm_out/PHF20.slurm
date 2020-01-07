#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PHF20.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PHF20 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PHF20L1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PHF21A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PHF21B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PHF23 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PHF24 
