#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CYB5RL.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CYB5RL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CYBA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CYBB 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CYBRD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CYC1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CYCS 
