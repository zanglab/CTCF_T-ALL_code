#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RAB41.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RAB41 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAB42 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAB43 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAB44 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAB4A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAB4B 
