#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o LGALS13.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i LGALS13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LGALS14 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LGALS2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LGALS3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LGALS3BP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LGALS4 
