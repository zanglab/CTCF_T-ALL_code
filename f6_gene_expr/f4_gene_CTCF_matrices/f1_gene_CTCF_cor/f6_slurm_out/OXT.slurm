#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o OXT.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i OXT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OXTR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i P2RX1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i P2RX2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i P2RX3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i P2RX4 
