#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ADORA2A.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ADORA2A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADORA2B 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADORA3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADPGK 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADPRH 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ADPRHL1 
