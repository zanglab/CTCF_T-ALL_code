#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SLC25A28.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC25A28 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC25A29 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC25A3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC25A30 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC25A31 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC25A32 
