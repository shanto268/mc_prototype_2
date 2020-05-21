c++  `root-config --cflags` -o muonAnalysis3 sc8muontree.cc muonAnalysis3_newsh.cc `root-config --glibs`
./muonAnalysis3  /Users/raulperezjr/hep/MuonSC8/mc/v4/sim/muonTree05_$1deg_$2simair.root  hist5_air$1deg.root
./muonAnalysis3  /Users/raulperezjr/hep/MuonSC8/mc/v4/sim/muonTree05_$1deg_$2simwater.root  hist5_water$1deg.root
