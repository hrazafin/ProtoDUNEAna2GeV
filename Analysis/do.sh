nEntryToStop=$1
which mkexe.sh

mkdir -p output

mkexe.sh Analysis 

./Analysis $nEntryToStop; 



