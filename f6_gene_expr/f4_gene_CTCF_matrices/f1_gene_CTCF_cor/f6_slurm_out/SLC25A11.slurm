#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SLC25A11.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC25A11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC25A12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC25A13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC25A14 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC25A15 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC25A16 
