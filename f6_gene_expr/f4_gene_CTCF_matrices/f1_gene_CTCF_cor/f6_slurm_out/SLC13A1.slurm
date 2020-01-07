#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SLC13A1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC13A1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC13A2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC13A3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC13A4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC13A5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLC14A1 
