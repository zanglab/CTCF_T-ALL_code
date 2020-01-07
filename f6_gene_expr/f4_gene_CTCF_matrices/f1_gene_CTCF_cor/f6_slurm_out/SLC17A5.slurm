#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SLC17A5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC17A5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC17A6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC17A7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC17A8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC17A9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC18A1 
