#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o UNC50.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i UNC50 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UNC5A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UNC5B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UNC5C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UNC5CL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i UNC5D 
