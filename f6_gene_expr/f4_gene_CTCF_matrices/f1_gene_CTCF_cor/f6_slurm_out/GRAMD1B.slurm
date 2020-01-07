#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GRAMD1B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GRAMD1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GRAMD1C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GRAMD4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GRAP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GRAP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GRAPL 
