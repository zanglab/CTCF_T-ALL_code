#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o C10orf67.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i C10orf67 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C10orf71 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C10orf76 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C10orf82 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C10orf88 
time python ../7_cor_CTCF_binding_GeneExpr.py -i C10orf90 
