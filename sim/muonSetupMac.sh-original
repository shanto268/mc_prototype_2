#!/bin/bash

#  setup ROOT environment
source /Applications/root_v6.18.00/bin/thisroot.sh

#  drop_from_path taken from /root/build/bin/thisroot.sh
drop_from_path()
{
   # Assert that we got enough arguments
   if test $# -ne 2 ; then
      echo "drop_from_path: needs 2 arguments"
      return 1
   fi

   local p=$1
   local drop=$2

   newpath=`echo $p | sed -e "s;:${drop}:;:;g" \
                          -e "s;:${drop}\$;;g"   \
                          -e "s;^${drop}:;;g"   \
                          -e "s;^${drop}\$;;g"`
}

export CRYHOME=/Users/raulperezjr/hep/cry/cry_v1.7
export CRYDATAPATH=/Users/raulperezjr/hep/cry/cry_v1.7/data

#  setup GEANT4 environment
. /Users/raulperezjr/hep/g4/Geant4-10.5.1-Darwin/bin/geant4.sh

export G4BASE=Users/raulperezjr/hep/g4/Geant4-10.5.1-Darwin/source
export G4INSTALL=/Users/raulperezjr/hep/g4/Geant4-10.5.1-Darwin/share/Geant4-10.5.1/geant4make

source /Users/raulperezjr/hep/g4/Geant4-10.5.1-Darwin/share/Geant4-10.5.1/geant4make/geant4make.sh
export G4BIN="$PWD"

export LD_LIBRARY_PATH=/Users/raulperezjr/hep/g4/Geant4-10.5.1-Darwin/lib
export DYLD_LIBRARY_PATH=/Users/raulperezjr/hep/g4/Geant4-10.5.1-Darwin/lib:/Applications/root_v6.18.00/lib
export SHLIB_PATH=/Users/raulperezjr/hep/g4/Geant4-10.5.1-Darwin/lib
export LIBPATH=/Users/raulperezjr/hep/g4/Geant4-10.5.1-Darwin/lib

