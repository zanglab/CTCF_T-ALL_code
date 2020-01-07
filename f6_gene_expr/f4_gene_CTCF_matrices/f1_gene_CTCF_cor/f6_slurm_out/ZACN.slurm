#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZACN.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZACN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZADH2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZAN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZAP70 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZAR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZAR1L 
