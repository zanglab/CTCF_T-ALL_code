#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZNF628.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF628 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF629 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF630 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF638 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF639 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF641 
