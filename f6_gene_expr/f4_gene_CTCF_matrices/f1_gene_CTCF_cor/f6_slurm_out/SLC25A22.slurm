#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SLC25A22.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC25A22 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC25A23 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC25A24 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC25A25 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC25A26 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC25A27 
