#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NES.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NES 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NET1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NETO1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NETO2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NEU1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NEU2 
