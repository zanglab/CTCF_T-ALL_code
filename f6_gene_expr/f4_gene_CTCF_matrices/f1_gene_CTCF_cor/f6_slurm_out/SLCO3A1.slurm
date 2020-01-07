#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SLCO3A1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SLCO3A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLCO4A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLCO4C1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLCO5A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLCO6A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLF1 
