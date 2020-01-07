#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o UCN2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i UCN2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UCN3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UCP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UCP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UCP3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UEVLD 
