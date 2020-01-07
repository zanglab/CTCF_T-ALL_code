#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SLC25A33.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC25A33 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC25A34 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC25A35 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC25A36 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC25A37 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC25A38 
