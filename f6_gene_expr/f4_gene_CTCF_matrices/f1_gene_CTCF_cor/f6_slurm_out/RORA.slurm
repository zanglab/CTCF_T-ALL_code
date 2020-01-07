#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RORA.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RORA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RORB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RORC 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ROS1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RP1L1 
