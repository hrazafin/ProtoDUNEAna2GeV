nEntryToStop=$1
which mkexe.sh

mkdir -p output

mkexe.sh Analysis -lRooUnfold -I${ROOUNFOLD}  -L${ROOUNFOLD} || return 1  

./Analysis $nEntryToStop; 



