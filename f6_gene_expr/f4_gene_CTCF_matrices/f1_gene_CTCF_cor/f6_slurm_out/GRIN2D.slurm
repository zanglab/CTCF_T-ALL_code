#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GRIN2D.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GRIN2D 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GRIN3A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GRIN3B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GRINA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GRIP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GRIP2 
