#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SLC46A3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC46A3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC47A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC47A2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC48A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC4A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC4A10 
