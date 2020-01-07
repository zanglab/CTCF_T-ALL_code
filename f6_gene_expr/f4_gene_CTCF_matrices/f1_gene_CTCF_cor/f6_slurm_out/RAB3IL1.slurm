#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RAB3IL1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RAB3IL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAB3IP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAB40A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAB40AL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAB40B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAB40C 
