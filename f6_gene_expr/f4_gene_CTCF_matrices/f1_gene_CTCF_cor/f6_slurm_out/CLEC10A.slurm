#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CLEC10A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CLEC10A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLEC11A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLEC12A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLEC12B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLEC14A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLEC16A 
