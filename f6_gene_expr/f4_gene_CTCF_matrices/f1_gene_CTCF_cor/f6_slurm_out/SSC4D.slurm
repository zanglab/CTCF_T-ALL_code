#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SSC4D.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SSC4D 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SSC5D 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SSFA2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SSH1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SSH2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SSH3 
