#!/bin/bash
#count the k-mers
mkdir -p  kmc_temp
./kmc -b -ci3 -fq -k32 -cs250 new.fastq list.kmc kmc_temp
./kmc_dump -ci3  list.kmc dump

#change the format of file called dump
python ../edit.py dump pair

#construct the read-colored matrix with reduced number of colors:
../reduce_color new.fastq 20000 32
