#!/bin/bash

tsp -S 8

cd ${GENIE}/config
cp Messenger_laconic.xml Messenger.xml 
cd -

DUNE_ND_FLUX=$(nuis flux DUNE_ND)
DUNE_ND_FLUX_RANGE=$(nuis flux range ${DUNE_ND_FLUX})

for i in {0..20}; do

  sleep 1

  #generate in separate subfolders to make sure they don't argue while
  mkdir gen_${i}; cd gen_${i}
  tsp gevgen -p 14 -t 1000180400 -n 1000000 \
    -f ${DUNE_ND_FLUX} -e ${DUNE_ND_FLUX_RANGE} \
    -o genie.DUNE_ND.numu.${i}.root \
    --tune ${GENIE_XSEC_TUNE} --cross-sections ${GENIE_XSEC_FILE}
  cd -

done

