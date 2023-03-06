which mkexe.sh

mkdir -p output

mkexe.sh NewXSstudy -lRooUnfold -I${ROOUNFOLD}  -L${ROOUNFOLD} || return 1

./NewXSstudy
