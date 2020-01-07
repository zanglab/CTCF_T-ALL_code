#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o IZUMO1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i IZUMO1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IZUMO1R 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IZUMO2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IZUMO3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i IZUMO4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i JADE1 
