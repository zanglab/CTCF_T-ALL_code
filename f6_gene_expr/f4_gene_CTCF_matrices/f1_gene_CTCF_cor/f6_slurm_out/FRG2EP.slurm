#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FRG2EP.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FRG2EP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FRK 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FRMD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FRMD3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FRMD4A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FRMD4B 
