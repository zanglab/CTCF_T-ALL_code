#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o MGAT5B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i MGAT5B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MGEA5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MGLL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MGME1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MGMT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i MGP 
