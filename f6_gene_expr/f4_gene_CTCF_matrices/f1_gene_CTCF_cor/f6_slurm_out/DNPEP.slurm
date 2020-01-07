#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o DNPEP.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i DNPEP 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DNPH1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DNTT 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DNTTIP1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DNTTIP2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i DOC2A 
