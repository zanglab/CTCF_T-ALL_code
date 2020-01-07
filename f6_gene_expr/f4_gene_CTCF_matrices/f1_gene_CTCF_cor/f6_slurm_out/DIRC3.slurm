#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DIRC3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DIRC3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DIS3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DIS3L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DIS3L2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DISC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DISP1 
