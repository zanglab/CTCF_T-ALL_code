#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GVINP1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GVINP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GXYLT1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GXYLT1P3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GXYLT2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GYG1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GYG2 
