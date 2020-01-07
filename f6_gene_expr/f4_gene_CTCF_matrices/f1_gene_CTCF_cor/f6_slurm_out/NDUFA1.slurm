#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o NDUFA1.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i NDUFA1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NDUFA10 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NDUFA11 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NDUFA12 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NDUFA13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i NDUFA2 
