#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KIAA1211L.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KIAA1211L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KIAA1217 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KIAA1257 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KIAA1324 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KIAA1324L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KIAA1328 
