#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o UTY.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i UTY 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UVRAG 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UVSSA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UXS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UXT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i VAC14 
