#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ATP1B4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ATP1B4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATP2A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATP2A2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATP2A3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATP2B1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATP2B2 
