#!/bin/bash

set -e

//c++  `root-config --cflags` -o anabin anaWTP_v01a.cc   `root-config --glibs`

//./anabin ../Data/v4frun_$1_muon_wtp.txt  histWTP_$1.root

//echo done!

//cp histWTP_$1.root root2tex/
//cd root2tex/

//echo In $PWD

//c++ `root-config --cflags` -o root2tex_v9 root2tex_v9.cc `root-config --glibs`

//histWTP_$1.root
hist6_50deg.root

pdflatex plots_histWTP_1.tex

open plots_histWTP_1.pdf

echo PDF created
cd ..
