#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o CLEC17A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i CLEC17A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLEC18A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLEC18B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLEC18C 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLEC19A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i CLEC1A 
