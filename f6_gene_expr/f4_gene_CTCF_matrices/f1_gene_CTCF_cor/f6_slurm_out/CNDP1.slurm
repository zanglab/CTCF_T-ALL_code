#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CNDP1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CNDP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CNDP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CNEP1R1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CNFN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CNGA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CNGA2 
