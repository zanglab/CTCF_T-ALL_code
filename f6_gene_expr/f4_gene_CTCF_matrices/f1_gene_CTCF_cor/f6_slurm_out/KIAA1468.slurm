#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KIAA1468.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KIAA1468 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KIAA1522 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KIAA1549 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KIAA1549L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KIAA1551 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KIAA1586 
