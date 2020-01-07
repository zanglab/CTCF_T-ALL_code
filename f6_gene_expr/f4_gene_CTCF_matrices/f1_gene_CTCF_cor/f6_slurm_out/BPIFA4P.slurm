#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o BPIFA4P.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i BPIFA4P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BPIFB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BPIFB2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BPIFB3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BPIFB4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BPIFB6 
