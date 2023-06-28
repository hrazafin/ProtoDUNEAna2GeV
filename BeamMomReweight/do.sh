nEntryToStop=$1
which mkexe.sh

mkdir -p output

mkexe.sh BeamMomReweight -lRooUnfold -I${ROOUNFOLD}  -L${ROOUNFOLD} || return 1  

./BeamMomReweight $nEntryToStop; 



