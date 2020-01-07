#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o EPHB2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i EPHB2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EPHB3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EPHB4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EPHB6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EPHX1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EPHX2 
