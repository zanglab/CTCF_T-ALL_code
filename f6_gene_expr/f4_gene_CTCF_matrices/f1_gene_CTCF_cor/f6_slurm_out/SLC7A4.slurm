#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SLC7A4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC7A4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC7A5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC7A5P1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC7A5P2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC7A6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC7A6OS 
