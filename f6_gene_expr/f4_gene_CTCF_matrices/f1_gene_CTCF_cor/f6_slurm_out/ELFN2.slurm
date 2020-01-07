#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ELFN2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ELFN2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ELK1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ELK2AP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ELK3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ELK4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ELL 
