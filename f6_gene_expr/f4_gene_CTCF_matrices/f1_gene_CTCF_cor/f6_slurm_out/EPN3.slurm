#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o EPN3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i EPN3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EPO 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EPOR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EPPIN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EPPIN-WFDC6 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EPPK1 
