#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PPM1J.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PPM1J 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPM1K 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPM1L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPM1M 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPM1N 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PPME1 
