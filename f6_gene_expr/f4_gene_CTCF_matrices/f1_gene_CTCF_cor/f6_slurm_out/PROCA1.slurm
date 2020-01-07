#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PROCA1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PROCA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PROCR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRODH 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PRODH2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PROK1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PROK2 
