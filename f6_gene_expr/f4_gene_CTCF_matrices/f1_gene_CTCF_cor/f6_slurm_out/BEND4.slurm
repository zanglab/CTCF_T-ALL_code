#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o BEND4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i BEND4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BEND5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BEND6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BEND7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BEST1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BEST2 
