#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FLI1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FLI1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FLII 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FLNA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FLNB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FLNC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FLOT1 
