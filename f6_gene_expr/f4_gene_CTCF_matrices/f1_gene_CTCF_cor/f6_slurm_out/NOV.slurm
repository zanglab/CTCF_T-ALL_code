#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NOV.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NOV 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NOVA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NOVA2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NOX1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NOX3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NOX4 
