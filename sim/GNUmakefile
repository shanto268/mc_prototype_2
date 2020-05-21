# $Id: GNUmakefile,v 1.2 2000-10-19 12:22:10 stanaka Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

$(warning This is sks makefile)

name := exampleB4a
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

.PHONY: all
all: lib bin

# ROOT support
CPPFLAGS += -I$(shell root-config --incdir --libs --cflags)
EXTRALIBS += $(shell root-config --glibs)
CPPFLAGS += -I/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/include/c++/v1/


$(info CPPFLAGS 2 is [${CPPFLAGS}])

# EXTRALIBS += -L$(CRYHOME)/lib -lCRY
CPPFLAGS  += -I$(CRYHOME)/src

include $(G4INSTALL)/config/binmake.gmk

$(info G4TARGET is [${G4TARGET}])
# $(info G4BIN is [${G4BIN}])
$(info G4INSTALL is [${G4INSTALL}])
# $(info QTHOME is [${QTHOME}])
# $(info G4LIB_BUILD_EXPAT is [${G4LIB_BUILD_EXPAT}])
# $(info G4SYSTEM is [${G4SYSTEM}])
# $(info G4LIB_BUILD_SHARED is [${G4LIB_BUILD_SHARED}])
# $(info G4LIB_BUILD_STATIC is [${G4LIB_BUILD_STATIC}])
# $(info CRYHOME is [${CRYHOME}])
$(info CPPFLAGS is [${CPPFLAGS}])
# $(info EXTRALIBS is [${EXTRALIBS}])
$(info CXX is [${CXX}])

visclean:
	rm -f g4*.prim g4*.eps g4*.wrl
	rm -f .DAWN_*

