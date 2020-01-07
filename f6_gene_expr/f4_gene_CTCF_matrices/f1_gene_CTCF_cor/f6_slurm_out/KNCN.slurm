#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KNCN.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KNCN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KNDC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KNG1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KNOP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KNSTRN 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KNTC1 
