#!/bin/bash

tsp -S 8

for i in 0 1 2 3 4 5 6 7; do

  sleep 1

  #generate in separate subfolders to make sure they don't argue while
  mkdir gen_${i}; cd gen_${i}
  tsp NUIS_QUIET=ON nuis-gen-GENIE -E DUNE_ND -n 200000 -t Ar -o genie.DUNE_ND.numu.${i}.root
  cd -

done

