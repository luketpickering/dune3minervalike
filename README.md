# dune3minervalike

Build the sample with `make` from within this directory with NUISANCE set up.

Generate events by running `./gen.sh`. This will launch 8 processes to generate events, check on their status with `tsp`. Once they have all finished, which should take about 20 minutes, (check with `tsp`) run `./gen_merge.sh` to merge them into a single file in the current directory, `genie.DUNE_ND.numu.merged.root`.

Can then run the sample like: `nuis comp genie.DUNE_ND.numu.merged.root -t GENIE -s DUNE_3D_MINERvALike -o DUNE_3D_MINERvALike.out.root -f`.

If you want to be able to modify zexp parameters, you can use: `nuiscomp -c throwcard.xml`

To make throws: `./throws.sh`, this also uses `tsp` and so you need to monitor `tsp` to make sure that all of the throws jobs have finished before running `throws_merge.sh` to merge them into a single output file with error bands. Check the bash source code in these simple scripts to see how they work.