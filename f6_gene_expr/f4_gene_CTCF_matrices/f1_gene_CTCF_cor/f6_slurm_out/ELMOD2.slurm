#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ELMOD2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ELMOD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ELMOD3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ELMSAN1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ELN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ELOF1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ELOVL1 
