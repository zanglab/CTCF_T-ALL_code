#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SLC12A4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC12A4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC12A5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC12A6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC12A7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC12A8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC12A9 
