#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RAB27A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RAB27A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAB27B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAB28 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAB29 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAB2A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAB2B 
