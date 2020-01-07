#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o OSGIN1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i OSGIN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OSGIN2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OSM 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OSMR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OSR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OSR2 
