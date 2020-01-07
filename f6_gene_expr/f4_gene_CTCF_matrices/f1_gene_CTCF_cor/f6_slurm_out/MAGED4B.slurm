#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MAGED4B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MAGED4B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAGEE1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAGEE2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAGEF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAGEH1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAGEL2 
