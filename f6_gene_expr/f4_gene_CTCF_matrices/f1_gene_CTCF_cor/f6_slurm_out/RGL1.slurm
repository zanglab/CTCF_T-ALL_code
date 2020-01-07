#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RGL1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RGL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RGL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RGL3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RGL4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RGMA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RGMB 
