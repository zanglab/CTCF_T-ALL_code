#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MRPS22.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MRPS22 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MRPS23 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MRPS24 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MRPS25 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MRPS26 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MRPS27 
