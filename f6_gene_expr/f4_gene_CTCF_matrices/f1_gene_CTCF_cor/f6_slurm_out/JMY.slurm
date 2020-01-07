#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o JMY.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i JMY 
time python ../7_cor_CTCF_binding_GeneExpr.py -i JOSD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i JOSD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i JPH1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i JPH2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i JPH3 
