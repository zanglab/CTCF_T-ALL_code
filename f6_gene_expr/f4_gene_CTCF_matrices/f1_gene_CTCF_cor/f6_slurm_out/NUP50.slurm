#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NUP50.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NUP50 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NUP54 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NUP58 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NUP62 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NUP62CL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NUP85 
