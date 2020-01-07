#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o LCE1E.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i LCE1E 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LCE1F 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LCE2A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LCE2B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LCE2C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LCE3A 
