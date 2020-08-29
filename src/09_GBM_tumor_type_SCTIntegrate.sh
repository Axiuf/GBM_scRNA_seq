#!/bin/bash

for integrate_dim in 10 20 30 50 70 100
do
  Rscript ./src/09_GBM_tumor_type_SCTIntegrate.R $integrate_dim
done

wait
