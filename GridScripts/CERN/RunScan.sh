#!/bin/bash

i=$1
j=$2
k=$3

source /cvmfs/sft.cern.ch/lcg/contrib/gcc/4.8/x86_64-slc6/setup.sh && source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.04.18/x86_64-slc6-gcc48-opt/root/bin/thisroot.sh
echo "ls before copying input files: "
ls --color -lh -c
cp /eos/user/j/jsotooto/DUNE/lightana/superclustering/newtree/ClusteringScan/macrogrid.C .
cp /eos/user/j/jsotooto/DUNE/lightana/superclustering/newtree/*.h .
cp /eos/user/j/jsotooto/DUNE/lightana/superclustering/newtree/*.C .
echo "ls after copying input files: "
ls --color -lh -c

root -l -b -q "macrogrid.C+g($i,$j,$k)"

echo "ls after running macro: "
ls --color -lh -c

cp *.root /eos/user/j/jsotooto/DUNE/lightana/data/HeavyRun/

echo "ls after deleting files: "
ls --color -lh -c


echo "done"
