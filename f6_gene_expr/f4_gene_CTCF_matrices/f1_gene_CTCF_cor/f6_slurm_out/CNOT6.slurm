#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CNOT6.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CNOT6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CNOT6L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CNOT7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CNOT8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CNP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CNPPD1 
