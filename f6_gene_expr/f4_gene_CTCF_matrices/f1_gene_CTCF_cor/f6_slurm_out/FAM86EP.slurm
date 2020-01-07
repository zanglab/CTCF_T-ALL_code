#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FAM86EP.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM86EP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM86FP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM86HP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM86JP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM89A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM89B 
