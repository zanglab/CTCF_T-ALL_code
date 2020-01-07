#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o HEPACAM.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i HEPACAM 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HEPACAM2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HEPH 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HEPHL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HEPN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HERC1 
