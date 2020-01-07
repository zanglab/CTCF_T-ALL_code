#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NRXN2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NRXN2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NRXN3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NSA2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NSD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NSDHL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NSF 
