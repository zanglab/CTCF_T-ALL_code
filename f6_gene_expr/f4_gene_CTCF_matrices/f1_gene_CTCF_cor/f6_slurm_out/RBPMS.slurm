#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RBPMS.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RBPMS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RBPMS2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RBSN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RBX1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RC3H1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RC3H2 
