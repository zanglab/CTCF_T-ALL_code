#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SLCO1B1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SLCO1B1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLCO1B3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLCO1B7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLCO1C1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLCO2A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLCO2B1 
