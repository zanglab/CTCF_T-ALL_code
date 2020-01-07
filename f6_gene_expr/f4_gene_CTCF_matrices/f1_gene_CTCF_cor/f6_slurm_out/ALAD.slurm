#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ALAD.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ALAD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ALAS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ALAS2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ALB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ALCAM 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ALDH16A1 
