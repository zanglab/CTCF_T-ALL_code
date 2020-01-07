#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RRS1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RRS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RSAD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RSAD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RSBN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RSBN1L 
