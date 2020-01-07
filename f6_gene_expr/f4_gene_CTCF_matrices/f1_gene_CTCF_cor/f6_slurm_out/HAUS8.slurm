#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o HAUS8.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i HAUS8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HAVCR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HAVCR1P1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HAVCR2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HAX1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HBA1 
