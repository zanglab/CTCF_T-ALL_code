#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MAGI1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MAGI1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAGI2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAGI3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAGIX 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAGOH 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAGOH2P 
