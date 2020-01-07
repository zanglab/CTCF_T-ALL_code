#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SLC39A6.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC39A6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC39A7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC39A8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC39A9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC3A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC3A2 
