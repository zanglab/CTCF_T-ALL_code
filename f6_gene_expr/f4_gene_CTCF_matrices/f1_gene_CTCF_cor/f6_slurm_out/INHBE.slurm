#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o INHBE.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i INHBE 
time python ../7_cor_CTCF_binding_GeneExpr.py -i INIP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i INMT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i INO80 
time python ../7_cor_CTCF_binding_GeneExpr.py -i INO80B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i INO80B-WBP1 
