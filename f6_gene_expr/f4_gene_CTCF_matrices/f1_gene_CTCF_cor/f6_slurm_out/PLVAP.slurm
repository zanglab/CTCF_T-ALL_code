#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PLVAP.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PLVAP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLXDC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLXDC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLXNA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLXNA2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PLXNA3 
