#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o WNT8A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i WNT8A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WNT8B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WNT9A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WNT9B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WRAP53 
time python ../7_cor_CTCF_binding_GeneExpr.py -i WRAP73 
