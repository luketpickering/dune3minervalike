# dune3minervalike

build with `make` from within this directory with NUISANCE set up.

Generate events by running `./gen.sh`. This will launch 8 processes to generate events, check on their status with `tsp`. Once they have all finished, which should take about 20 minutes, (check with `tsp`) run `./gen_merge.sh` to merge them into a single file in the current directory, `genie.DUNE_ND.numu.merged.root`.

Can then run the sample like: `nuis comp genie.DUNE_ND.numu.merged.root -t GENIE -s DUNE_3D_MINERvALike -o DUNE_3D_MINERvALike.out.root -f`.