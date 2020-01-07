#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SLK.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SLK 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLMAP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLPI 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLTM 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLU7 
