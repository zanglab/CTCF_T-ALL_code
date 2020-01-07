#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o LRRFIP1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i LRRFIP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRRFIP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRRIQ1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRRIQ3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRRIQ4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LRRK1 
