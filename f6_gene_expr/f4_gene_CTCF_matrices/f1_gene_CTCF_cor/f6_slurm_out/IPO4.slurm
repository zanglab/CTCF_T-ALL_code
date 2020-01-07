#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o IPO4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i IPO4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IPO5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IPO5P1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IPO7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IPO8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IPO9 
