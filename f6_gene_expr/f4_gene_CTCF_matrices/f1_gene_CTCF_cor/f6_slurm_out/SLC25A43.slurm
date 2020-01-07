#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SLC25A43.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC25A43 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC25A44 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC25A45 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC25A46 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC25A47 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC25A48 
