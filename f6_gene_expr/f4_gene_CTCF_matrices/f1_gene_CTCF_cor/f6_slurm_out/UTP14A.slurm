#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o UTP14A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i UTP14A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UTP14C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UTP15 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UTP18 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UTP20 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UTP23 
