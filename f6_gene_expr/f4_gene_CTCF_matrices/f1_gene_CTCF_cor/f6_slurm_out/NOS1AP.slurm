#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NOS1AP.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NOS1AP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NOS2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NOS3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NOSIP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NOSTRIN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NOTCH1 
