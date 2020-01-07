#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o POLE2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i POLE2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POLE3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POLE4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POLG 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POLG2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i POLH 
