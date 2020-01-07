#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o P2RX5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i P2RX5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i P2RX5-TAX1BP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i P2RX6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i P2RX6P 
time python ../7_cor_CTCF_binding_GeneExpr.py -i P2RX7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i P2RY1 
