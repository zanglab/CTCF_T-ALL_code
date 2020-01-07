#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NPEPPS.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NPEPPS 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NPFF 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NPFFR1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NPFFR2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NPHP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NPHP3 
