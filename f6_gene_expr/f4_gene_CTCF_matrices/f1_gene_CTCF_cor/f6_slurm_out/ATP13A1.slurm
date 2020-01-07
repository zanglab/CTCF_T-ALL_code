#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ATP13A1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ATP13A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATP13A2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATP13A3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATP13A4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATP13A5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATP1A1 
