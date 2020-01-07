#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NAIF1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NAIF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NAIP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NALCN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NAMPT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NANOG 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NANOGNB 
