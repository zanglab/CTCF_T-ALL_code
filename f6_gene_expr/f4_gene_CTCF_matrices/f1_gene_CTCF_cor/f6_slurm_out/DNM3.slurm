#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DNM3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DNM3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DNMBP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DNMT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DNMT3A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DNMT3B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DNMT3L 
