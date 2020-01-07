#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o URM1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i URM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UROC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UROD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UROS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i USB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i USE1 
