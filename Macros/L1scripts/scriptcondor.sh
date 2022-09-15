#!/bin/bash

cd ~/Run3DQM/CMSSW_12_4_3/src/JetMETStudies/Macros/L1scripts/
cmsenv 


python3 analyzel1.py --maxEvents -1 -i $1 -o $2 -c $3

