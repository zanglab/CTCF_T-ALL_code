#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RDH13.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RDH13 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RDH14 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RDH16 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RDH5 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RDH8 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RDM1 
