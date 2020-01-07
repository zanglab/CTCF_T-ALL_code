#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NUGGC.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NUGGC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NUMA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NUMB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NUMBL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NUP107 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NUP133 
