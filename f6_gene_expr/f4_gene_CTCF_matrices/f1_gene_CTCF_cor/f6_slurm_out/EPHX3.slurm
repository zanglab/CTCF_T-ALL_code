#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o EPHX3.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i EPHX3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EPHX4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EPM2A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EPM2AIP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EPN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i EPN2 
