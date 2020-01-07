#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FADS6.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FADS6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAF2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAH 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAHD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAHD2A 
