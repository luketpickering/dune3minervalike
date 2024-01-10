#!/bin/bash

NTHROWS=150

for i in 0 1 2 3 4 5 6 7; do

  sleep 1

  mkdir gen_${i}; cd gen_${i}
  tsp nuissyst -c throwcard.xml -f ThrowErrors -o dunesyst.${i} \
    -q error_throws=${NTHROWS} -s $(( NTHROWS * i )) -n 500000
  cd -

done