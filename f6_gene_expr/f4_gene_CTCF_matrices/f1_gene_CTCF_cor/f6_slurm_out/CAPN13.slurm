#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CAPN13.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CAPN13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CAPN14 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CAPN15 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CAPN2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CAPN3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CAPN5 
