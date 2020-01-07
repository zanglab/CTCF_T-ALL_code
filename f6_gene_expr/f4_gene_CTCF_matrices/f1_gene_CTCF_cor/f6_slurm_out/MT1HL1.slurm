#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MT1HL1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MT1HL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MT1JP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MT1L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MT1M 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MT1X 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MT2A 
