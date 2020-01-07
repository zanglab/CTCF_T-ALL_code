#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CDK14.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CDK14 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDK15 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDK16 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDK17 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDK18 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDK19 
