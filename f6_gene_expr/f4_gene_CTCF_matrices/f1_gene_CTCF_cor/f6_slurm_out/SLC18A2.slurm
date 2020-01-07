#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SLC18A2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC18A2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC18A3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC18B1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC19A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC19A2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC19A3 
