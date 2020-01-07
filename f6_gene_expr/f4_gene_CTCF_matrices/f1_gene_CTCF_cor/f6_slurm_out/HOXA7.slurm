#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o HOXA7.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i HOXA7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HOXA9 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HOXB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HOXB13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HOXB2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HOXB3 
