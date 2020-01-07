#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o LMX1B.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i LMX1B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LNP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LNPEP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LNX1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LNX2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LONP1 
