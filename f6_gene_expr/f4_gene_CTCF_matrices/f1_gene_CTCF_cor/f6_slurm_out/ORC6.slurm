#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ORC6.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ORC6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ORM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ORM2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ORMDL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ORMDL2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ORMDL3 
