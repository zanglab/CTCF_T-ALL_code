#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o FAM219B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM219B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM21EP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM220A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM221A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM221B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i FAM222A 
