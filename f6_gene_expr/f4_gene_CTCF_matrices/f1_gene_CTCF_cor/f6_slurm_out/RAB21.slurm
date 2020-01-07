#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o RAB21.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i RAB21 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAB22A 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAB23 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAB24 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAB25 
time python ../7_cor_CTCF_binding_GeneExpr.py -i RAB26 
