#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o UNC13A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i UNC13A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UNC13B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UNC13C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UNC13D 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UNC45A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UNC45B 
