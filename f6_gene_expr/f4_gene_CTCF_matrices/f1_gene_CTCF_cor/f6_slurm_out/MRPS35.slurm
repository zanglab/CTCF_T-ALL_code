#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MRPS35.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MRPS35 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MRPS36 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MRPS5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MRPS6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MRPS7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MRPS9 
