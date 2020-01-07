#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CARF.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CARF 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CARHSP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CARM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CARNMT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CARNS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CARS 
