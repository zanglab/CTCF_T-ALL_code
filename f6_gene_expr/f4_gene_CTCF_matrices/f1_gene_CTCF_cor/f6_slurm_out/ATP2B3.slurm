#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ATP2B3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ATP2B3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATP2B4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATP2C1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATP2C2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATP4A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATP4B 
