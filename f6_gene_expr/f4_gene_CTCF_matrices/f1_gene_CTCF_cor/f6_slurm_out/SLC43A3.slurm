#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SLC43A3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC43A3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC44A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC44A2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC44A3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC44A4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC44A5 
