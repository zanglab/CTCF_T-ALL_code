#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o LACRT.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i LACRT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LACTB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LACTB2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LACTBL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LAD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LAG3 
