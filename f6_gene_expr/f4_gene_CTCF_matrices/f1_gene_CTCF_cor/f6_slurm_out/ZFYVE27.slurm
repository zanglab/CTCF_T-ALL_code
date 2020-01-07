#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZFYVE27.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZFYVE27 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZFYVE28 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZFYVE9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZG16 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZG16B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZGLP1 
