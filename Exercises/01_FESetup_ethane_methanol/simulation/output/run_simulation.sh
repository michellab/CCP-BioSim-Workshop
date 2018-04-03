#!/bin/bash

for l in 0.00 0.30 0.60 1.00
do
cd lambda-$l
somd-freenrg -C ../../input/sim.cfg -c ../../input/vacuum.rst7 -t ../../input/vacuum.parm7 -m ../../input/MORPH.pert -p CPU -l $l
cd ../
done 
