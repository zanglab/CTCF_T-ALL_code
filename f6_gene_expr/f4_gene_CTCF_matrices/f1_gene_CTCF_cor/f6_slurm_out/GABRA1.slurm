#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o GABRA1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i GABRA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GABRA2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GABRA3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GABRA4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GABRA5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i GABRA6 
