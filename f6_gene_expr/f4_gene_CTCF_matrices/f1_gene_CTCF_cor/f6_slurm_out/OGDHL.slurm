#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o OGDHL.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i OGDHL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OGFOD1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OGFOD2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OGFOD3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OGFR 
time python ../7_cor_CTCF_binding_GeneExpr.py -i OGFRL1 
