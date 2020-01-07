#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CDH17.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CDH17 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDH18 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDH19 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDH2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDH20 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDH22 
