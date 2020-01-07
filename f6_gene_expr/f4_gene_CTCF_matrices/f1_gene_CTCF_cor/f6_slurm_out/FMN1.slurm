#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FMN1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FMN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FMN2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FMNL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FMNL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FMNL3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FMO1 
