#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o HEBP2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i HEBP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HECA 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HECTD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HECTD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HECTD3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i HECTD4 
