#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DPEP1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DPEP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DPEP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DPEP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DPF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DPF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DPF3 
