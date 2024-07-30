#!/bin/bash

# run in the data folder for this project
# on biowulf2
mkdir -p 00log

sbcmd="sbatch --cpus-per-task={threads} \
--mem={cluster.mem} \
--time={cluster.time} \
--partition={cluster.partition} \
--output={cluster.output} \
--error={cluster.error} \
{cluster.extra}"


snakemake -s $1 \
  --rerun-triggers mtime \
  -pr --local-cores 2 --jobs 500 \
  --cluster-config $3 \
  --cluster "$sbcmd"  --latency-wait 120 --rerun-incomplete \
  --configfile $2 --use-conda \
  -k --restart-times 0 
