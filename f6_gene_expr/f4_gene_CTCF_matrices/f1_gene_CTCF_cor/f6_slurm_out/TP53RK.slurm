#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TP53RK.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TP53RK 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TP53TG3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TP53TG3B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TP53TG3C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TP53TG3D 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TP53TG5 
