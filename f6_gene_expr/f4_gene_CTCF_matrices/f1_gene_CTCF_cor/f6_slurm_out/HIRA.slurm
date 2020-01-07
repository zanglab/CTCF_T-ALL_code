#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o HIRA.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i HIRA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HIRIP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HIST1H1A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HIST1H1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HIST1H1C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HIST1H1D 
