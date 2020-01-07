#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MARVELD2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MARVELD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MARVELD3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MAS1L 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MASP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MASP2 
