#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o BACE1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i BACE1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BACE2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BACH1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BACH2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BAD 
time python ../7_cor_CTCF_binding_GeneExpr.py -i BAG1 
