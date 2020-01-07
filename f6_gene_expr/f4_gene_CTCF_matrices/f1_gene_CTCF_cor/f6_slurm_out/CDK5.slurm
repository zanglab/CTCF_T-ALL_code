#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CDK5.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CDK5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDK5R1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDK5R2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDK5RAP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDK5RAP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CDK5RAP3 
