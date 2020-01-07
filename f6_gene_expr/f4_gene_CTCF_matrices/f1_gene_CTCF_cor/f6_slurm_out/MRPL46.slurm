#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MRPL46.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MRPL46 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MRPL47 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MRPL48 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MRPL49 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MRPL50 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MRPL51 
