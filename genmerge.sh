# #!/bin/bash

# FILES="/home/abipeake/scratch/Formfactorstuff/Feb24022/"
# for i in {0..24}; do
#   FILES="${FILES},gen_${i}/genie.DUNE_ND.numu.${i}.root"
  
# done

# PrepareGENIE -i ${FILES} -f $(nuis flux DUNE_ND) -t 1000180400[1] -o genie.DUNE_ND.numu.merged.root


#!/bin/bash

DIR="/scratch/abipeake/Formfactorstuff/Feb24022"
FILES=""
for i in {0..24}; do
  FILES="${FILES},${DIR}/gen_${i}/genie.DUNE_ND.numu.${i}.root"
done
FILES="${FILES#,}"    # remove the leading comma

PrepareGENIE -i ${FILES} -f $(nuis flux DUNE_ND) -t "1000180400[1]" -o genie.DUNE_ND.numu.merged.root