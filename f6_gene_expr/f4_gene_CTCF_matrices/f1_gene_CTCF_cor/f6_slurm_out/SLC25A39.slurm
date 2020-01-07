#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SLC25A39.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC25A39 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC25A3P1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC25A4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC25A40 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC25A41 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC25A42 
