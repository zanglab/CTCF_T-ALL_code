#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RASD2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RASD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RASEF 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RASGEF1A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RASGEF1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RASGEF1C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RASGRF1 
