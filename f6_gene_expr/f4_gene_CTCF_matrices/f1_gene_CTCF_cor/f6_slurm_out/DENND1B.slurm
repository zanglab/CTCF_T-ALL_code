#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DENND1B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DENND1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DENND1C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DENND2A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DENND2C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DENND2D 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DENND3 
