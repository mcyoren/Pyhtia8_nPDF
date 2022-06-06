# PYTHIA configuration file.
# Generated on Пт 14 янв 2022 19:10:56 MSK with the user supplied options:
# --with-root=/usr/root
# --with-lhapdf6=/home/yoren/bnl/pythia8303/LHAPDF-6.4.0/install

# Install directory prefixes.
PREFIX_BIN=/home/yoren/bnl/pythia8303/bin
PREFIX_INCLUDE=/home/yoren/bnl/pythia8303/include
PREFIX_LIB=/home/yoren/bnl/pythia8303/lib
PREFIX_SHARE=/home/yoren/bnl/pythia8303/share/Pythia8

# Compilation flags (see ./configure --help for further documentation).
CXX=g++
CXX_COMMON=-O2 -std=c++11 -pedantic -W -Wall -Wshadow -fPIC
CXX_SHARED=-shared
CXX_SONAME=-Wl,-soname,
LIB_SUFFIX=.so
OBJ_COMMON=

EVTGEN_USE=false
EVTGEN_CONFIG=
EVTGEN_BIN=
EVTGEN_INCLUDE=
EVTGEN_LIB=

FASTJET3_USE=false
FASTJET3_CONFIG=
FASTJET3_BIN=
FASTJET3_INCLUDE=
FASTJET3_LIB=

HEPMC2_USE=false
HEPMC2_CONFIG=
HEPMC2_BIN=
HEPMC2_INCLUDE=
HEPMC2_LIB=

HEPMC3_USE=false
HEPMC3_CONFIG=
HEPMC3_BIN=
HEPMC3_INCLUDE=
HEPMC3_LIB=

LHAPDF5_USE=false
LHAPDF5_CONFIG=
LHAPDF5_BIN=
LHAPDF5_INCLUDE=
LHAPDF5_LIB=

LHAPDF6_USE=true
LHAPDF6_CONFIG=lhapdf-config
LHAPDF6_BIN=/home/yoren/bnl/pythia8303/LHAPDF-6.4.0/install/bin/
LHAPDF6_INCLUDE=-I/home/yoren/bnl/pythia8303/LHAPDF-6.4.0/install/include
LHAPDF6_LIB=-L/home/yoren/bnl/pythia8303/LHAPDF-6.4.0/install/lib -Wl,-rpath,/home/yoren/bnl/pythia8303/LHAPDF-6.4.0/install/lib -lLHAPDF

POWHEG_USE=false
POWHEG_CONFIG=lhapdf-config
POWHEG_BIN=
POWHEG_INCLUDE=
POWHEG_LIB=

RIVET_USE=false
RIVET_CONFIG=lhapdf-config
RIVET_BIN=
RIVET_INCLUDE=
RIVET_LIB=

ROOT_USE=true
ROOT_CONFIG=root-config
ROOT_BIN=/usr/root/bin/
ROOT_INCLUDE=-I/usr/root/include
ROOT_LIB=-L/usr/root/lib -Wl,-rpath,/usr/root/lib -lCore

YODA_USE=false
YODA_CONFIG=root-config
YODA_BIN=
YODA_INCLUDE=
YODA_LIB=

GZIP_USE=false
GZIP_CONFIG=root-config
GZIP_BIN=
GZIP_INCLUDE=
GZIP_LIB=

PYTHON_USE=false
PYTHON_CONFIG=root-config
PYTHON_BIN=
PYTHON_INCLUDE=
PYTHON_LIB=

MG5MES_USE=false
MG5MES_CONFIG=root-config
MG5MES_BIN=
MG5MES_INCLUDE=
MG5MES_LIB=

OPENMP_USE=false
OPENMP_CONFIG=root-config
OPENMP_BIN=
OPENMP_INCLUDE=
OPENMP_LIB=
