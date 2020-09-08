#!/bin/bash

#  setup ROOT environment

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

export CRYHOME=~/hep/cry/cry_v1.7
export CRYDATAPATH=~/hep/cry/cry_v1.7/data

export G4BASE=~/hep/g4/geant4.10.06.p01-install/source
export G4INSTALL=~/hep/g4/geant4.10.06.p01-install/share/Geant4-10.6.1/geant4make

source ~/hep/g4/geant4.10.06.p01-install/share/Geant4-10.6.1/geant4make/geant4make.sh
export G4BIN="$PWD"

