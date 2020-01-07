#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GALNT10.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GALNT10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GALNT11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GALNT12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GALNT13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GALNT14 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GALNT15 
