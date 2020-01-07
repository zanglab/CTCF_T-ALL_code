#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MBTPS1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MBTPS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MBTPS2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MC1R 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MC2R 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MC4R 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MC5R 
