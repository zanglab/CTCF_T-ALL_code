#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o IRS4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i IRS4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IRX1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IRX2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IRX3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IRX4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IRX5 
