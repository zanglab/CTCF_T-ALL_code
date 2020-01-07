#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o AHCYL1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i AHCYL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AHCYL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AHDC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AHI1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AHNAK 
time python ../7_cor_CTCF_binding_GeneExpr.py -i AHNAK2 
