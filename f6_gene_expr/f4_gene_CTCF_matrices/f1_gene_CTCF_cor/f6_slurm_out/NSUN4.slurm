#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NSUN4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NSUN4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NSUN5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NSUN5P1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NSUN5P2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NSUN6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NSUN7 
