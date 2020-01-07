#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o BMP2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i BMP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BMP2K 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BMP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BMP4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BMP5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BMP6 
