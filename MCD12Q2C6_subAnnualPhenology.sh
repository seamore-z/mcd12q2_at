#!/bin/bash
echo Submitting tile $1
echo Processing year $2
R --vanilla < /projectnb/modislc/users/seamorez/MCD12Q2/MCD12Q2C6_AnnualPhenology.R --args -tile $1 -year $2

