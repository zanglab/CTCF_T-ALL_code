#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o TMTC2.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i TMTC2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMTC3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMTC4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMUB1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMUB2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i TMX1 
