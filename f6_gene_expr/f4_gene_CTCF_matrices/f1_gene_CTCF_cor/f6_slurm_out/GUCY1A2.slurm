#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GUCY1A2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GUCY1A2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GUCY1B2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GUCY2C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GUCY2D 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GUCY2EP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GUCY2F 
