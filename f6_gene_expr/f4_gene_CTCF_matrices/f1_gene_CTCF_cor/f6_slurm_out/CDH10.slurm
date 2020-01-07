#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CDH10.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CDH10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDH11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDH12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDH13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDH15 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDH16 
