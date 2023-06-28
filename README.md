# How to get started

Find your favourite base directory and clone this repository
```
git clone https://github.com/GitKangaroo/protoDUNEAnalysis.git
```
Set up the environment 
```
cd protoDUNEAnalysis
source setup.sh
```
Also need the RooUnfold package to make the codes working
```
git clone https://github.com/roofit-dev/RooUnfold.git
cd RooUnfold
mkdir build
cd build
cmake ..
make -j4
```
Now, go to the main analysis folder and link the input files (dunegpvm server)
```
cd ../../Analysis/input
rm protoDUNE_mc_reco_flattree_prod4a_ntuple.root
ln -s /dune/data/users/kyang/pduneana_production/protoDUNE_mc_reco_flattree_prod4a_ntuple.root
rm protoDUNE_data_reco_flattree_prod4_ntuple.root
ln -s /dune/data/users/kyang/pduneana_production/protoDUNE_data_reco_flattree_prod4_ntuple.root
```
You are free to go, just type to run the full data and MC samples (takes ~ 15 mins)
```
./do.sh > log.txt
```
You can also specify the number of entries for a small sample just for a quick check
```
./do.sh 30000 > log_test.txt
```
