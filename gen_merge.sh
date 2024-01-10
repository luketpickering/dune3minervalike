#!/bin/bash

FILES=""
for i in 0 1 2 3 4 5 6 7; do
  FILES="${FILES},gen_${i}/genie.DUNE_ND.numu.${i}.root"
done

PrepareGENIE -i ${FILES} -f $(nuis flux DUNE_ND) -t 1000180400[1] -o genie.DUNE_ND.numu.merged.root