#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CDH23.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CDH23 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDH24 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDH26 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDH3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDH4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDH5 
