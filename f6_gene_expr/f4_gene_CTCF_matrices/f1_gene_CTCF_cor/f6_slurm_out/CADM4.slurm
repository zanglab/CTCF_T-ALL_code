#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CADM4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CADM4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CADPS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CADPS2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CAGE1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CALB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CALB2 
