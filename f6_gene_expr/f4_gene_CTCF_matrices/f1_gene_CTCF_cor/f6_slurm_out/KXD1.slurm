#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KXD1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KXD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KY 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KYNU 
time python ../7_cor_CTCF_binding_GeneExpr.py -i L1CAM 
time python ../7_cor_CTCF_binding_GeneExpr.py -i L1TD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i L2HGDH 
