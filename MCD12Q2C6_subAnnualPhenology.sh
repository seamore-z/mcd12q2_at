#!/bin/bash
echo Submitting tile $1
echo Processing year $2
R --vanilla < /projectnb/modislc/users/seamorez/mcd12q2_at/MCD12Q2C6_AnnualPhenology.R --args -tile $1 -year $2

