#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o OBSL1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i OBSL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OC90 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OCA2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OCEL1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OCIAD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OCIAD2 
