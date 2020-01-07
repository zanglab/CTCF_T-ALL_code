#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ADGRE3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ADGRE3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADGRE4P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADGRE5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADGRF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADGRF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADGRF3 
