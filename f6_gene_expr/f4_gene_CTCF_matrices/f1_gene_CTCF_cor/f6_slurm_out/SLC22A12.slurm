#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SLC22A12.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC22A12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC22A13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC22A14 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC22A15 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC22A16 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC22A17 
