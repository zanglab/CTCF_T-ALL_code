#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PIGC.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PIGC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PIGF 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PIGG 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PIGH 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PIGK 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PIGL 
