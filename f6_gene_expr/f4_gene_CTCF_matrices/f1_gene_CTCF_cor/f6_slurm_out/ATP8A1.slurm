#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ATP8A1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ATP8A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATP8A2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATP8B1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATP8B2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATP8B3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ATP8B4 
