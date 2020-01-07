#!/bin/bash

#SBATCH -n 1
#SBATCH --mem=4000
#SBATCH -t 24:00:00
#SBATCH -p standard
#SBATCH -A zanglab
#SBATCH -o LGALSL.out

time python ../7_cor_CTCF_binding_GeneExpr.py -i LGALSL 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LGI1 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LGI2 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LGI3 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LGI4 
time python ../7_cor_CTCF_binding_GeneExpr.py -i LGMN 
