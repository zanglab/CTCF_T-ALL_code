#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CDK1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CDK1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDK10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDK11A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDK11B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDK12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDK13 
