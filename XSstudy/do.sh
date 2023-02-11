which mkexe.sh

mkdir -p output

mkexe.sh XSstudy -lRooUnfold -I${ROOUNFOLD}  -L${ROOUNFOLD} || return 1

./XSstudy
