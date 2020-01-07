#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SH2B3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SH2B3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SH2D1A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SH2D1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SH2D2A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SH2D3A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SH2D3C 
