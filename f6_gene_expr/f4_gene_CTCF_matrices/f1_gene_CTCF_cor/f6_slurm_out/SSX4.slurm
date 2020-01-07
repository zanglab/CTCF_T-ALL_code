#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SSX4.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SSX4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SSX5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SSX6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SSX7 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SSX8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SSX9 
