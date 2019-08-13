#!/bin/bash

### Setup of CMSSW release
NUM_CORES=10
CMSSW=CMSSW_10_2_16_UL

export SCRAM_ARCH=slc7_amd64_gcc700
export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
source $VO_CMS_SW_DIR/cmsset_default.sh

scramv1 project $CMSSW; pushd $CMSSW/src
eval `scramv1 runtime -sh`

# combine on 102X slc7
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
cd HiggsAnalysis/CombinedLimit
git fetch origin
git checkout v8.0.1
cd -

# CombineHarvester (current fork of official master for 102X with fixes for BinByBin)
git clone https://github.com/KIT-CMS/CombineHarvester.git CombineHarvester

# SM analysis specific code
git clone https://github.com/KIT-CMS/SMRun2Legacy.git CombineHarvester/SMRun2Legacy

# compile everything
scramv1 b clean; scramv1 b -j $NUM_CORES
