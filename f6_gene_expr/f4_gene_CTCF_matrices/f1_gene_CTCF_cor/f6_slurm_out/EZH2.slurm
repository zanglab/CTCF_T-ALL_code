#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o EZH2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i EZH2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EZR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i F10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i F11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i F11R 
time python ../7_cor_CTCF_binding_GeneExpr.py -i F12 
