#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SEMA3C.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SEMA3C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SEMA3D 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SEMA3E 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SEMA3F 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SEMA3G 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SEMA4A 
