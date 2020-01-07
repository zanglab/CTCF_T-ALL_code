#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MYH10.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MYH10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MYH11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MYH13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MYH14 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MYH15 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MYH16 
