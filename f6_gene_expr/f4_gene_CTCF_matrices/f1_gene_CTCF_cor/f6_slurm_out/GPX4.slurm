#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GPX4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GPX4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPX5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPX6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPX7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPX8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GRAMD1A 
