#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CD226.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CD226 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CD24 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CD244 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CD247 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CD248 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CD27 
