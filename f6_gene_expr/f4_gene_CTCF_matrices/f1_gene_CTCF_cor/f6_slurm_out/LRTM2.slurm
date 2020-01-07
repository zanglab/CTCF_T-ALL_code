#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o LRTM2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i LRTM2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRTOMT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRWD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LSAMP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LSG1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LSM1 
