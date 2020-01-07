#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o METTL11B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i METTL11B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i METTL13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i METTL14 
time python ../7_cor_CTCF_binding_GeneExpr.py -i METTL15 
time python ../7_cor_CTCF_binding_GeneExpr.py -i METTL16 
time python ../7_cor_CTCF_binding_GeneExpr.py -i METTL17 
