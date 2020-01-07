#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KIAA0825.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KIAA0825 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KIAA0895 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KIAA0895L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KIAA0930 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KIAA1024 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KIAA1024L 
