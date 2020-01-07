#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SLC25A17.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC25A17 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC25A18 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC25A19 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC25A2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC25A20 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC25A21 
