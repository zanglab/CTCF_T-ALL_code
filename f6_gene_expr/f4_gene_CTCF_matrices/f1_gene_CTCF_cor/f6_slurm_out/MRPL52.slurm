#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MRPL52.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MRPL52 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MRPL53 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MRPL54 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MRPL55 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MRPL57 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MRPL9 
