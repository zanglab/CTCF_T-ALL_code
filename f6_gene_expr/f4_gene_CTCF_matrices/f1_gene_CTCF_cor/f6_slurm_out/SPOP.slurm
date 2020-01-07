#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SPOP.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SPOP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPOPL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPPL2A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SPPL2B 
