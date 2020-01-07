#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SLC22A3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC22A3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC22A31 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC22A4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC22A5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC22A6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC22A7 
