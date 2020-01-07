#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o LCE3D.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i LCE3D 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LCE3E 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LCE5A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LCE6A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LCK 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LCLAT1 
