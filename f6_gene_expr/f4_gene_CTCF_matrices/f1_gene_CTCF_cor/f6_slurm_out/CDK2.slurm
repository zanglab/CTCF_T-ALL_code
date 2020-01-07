#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CDK2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CDK2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDK20 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDK2AP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDK2AP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDK3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDK4 
