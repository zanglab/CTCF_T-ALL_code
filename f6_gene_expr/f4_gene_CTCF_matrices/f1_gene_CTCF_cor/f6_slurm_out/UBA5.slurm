#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o UBA5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i UBA5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UBA52 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UBA6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UBA7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UBAC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UBAC2 
