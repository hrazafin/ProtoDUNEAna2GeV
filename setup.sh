echo calling DUNEprotoTKI/setup.sh 

source DUNE.sh
export DUNEPROTOTKI=$(pwd)
export LOCALBIN=${DUNEPROTOTKI}/bin
export ROOUNFOLD=${DUNEPROTOTKI}/RooUnfold/build
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${ROOUNFOLD}
export PATH=$PATH:${LOCALBIN}
setup_fnal_security
#source DUNE.sh
