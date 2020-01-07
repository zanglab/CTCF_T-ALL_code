#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ADAM21.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ADAM21 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADAM21P1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADAM22 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADAM23 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADAM28 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADAM29 
