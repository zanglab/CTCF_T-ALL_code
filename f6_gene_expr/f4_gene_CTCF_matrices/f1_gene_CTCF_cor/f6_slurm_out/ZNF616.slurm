#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o ZNF616.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF616 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF618 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF619 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF620 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF621 
time python ../7_cor_CTCF_binding_GeneExpr.py -i ZNF622 
