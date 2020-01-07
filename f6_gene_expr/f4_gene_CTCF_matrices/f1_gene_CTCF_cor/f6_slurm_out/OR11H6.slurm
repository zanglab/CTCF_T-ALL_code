#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o OR11H6.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i OR11H6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR11L1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR12D2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR12D3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR13A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OR13C2 
