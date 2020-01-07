#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o KCTD16.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i KCTD16 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCTD17 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCTD18 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCTD19 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCTD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i KCTD20 
