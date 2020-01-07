#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o PEAR1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i PEAR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PEBP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PEBP4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PECAM1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PECR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i PEF1 
