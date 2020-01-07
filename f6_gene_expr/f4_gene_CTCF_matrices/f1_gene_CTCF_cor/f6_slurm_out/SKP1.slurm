#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o SKP1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i SKP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SKP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLA2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLAIN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i SLAIN2 
