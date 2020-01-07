#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MEPE.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MEPE 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MERTK 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MESP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MESP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MEST 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MET 
