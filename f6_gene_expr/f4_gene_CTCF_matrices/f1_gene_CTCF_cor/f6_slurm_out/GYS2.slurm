#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GYS2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GYS2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GZF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GZMA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GZMB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GZMH 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GZMK 
