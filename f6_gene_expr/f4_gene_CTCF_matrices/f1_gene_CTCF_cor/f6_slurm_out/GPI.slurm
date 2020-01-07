#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GPI.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GPI 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPIHBP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPKOW 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPLD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPM6A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GPM6B 
